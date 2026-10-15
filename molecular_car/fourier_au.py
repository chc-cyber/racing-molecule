"""Conservative, bounded Fourier proxy of the frozen six-layer SAIP surface.

Only substrate queries are periodic. Direction dependence is retained through
the chart a=nx/nz, b=ny/nz and tensor bicubic Hermite interpolation. The chart
is a finite approximation, not a quadratic normal model or Bloch wavefunction.
Outside its support the potential smoothly returns to the validated lattice
sum, including derivatives of every switching weight.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
import zipfile

import numpy as np

from .periodic_au import (PeriodicAuPotential, ROOT, _taper, _anisotropic,
    _normal_geometry, _add_normal_chain_forces, _cache_weight, _radial,
    _cubic_weights, NORMAL_AREA_EPSILON, BONDED_NORMAL_POLICY)

ALGORITHM_VERSION = "fourier-saip-hermite-1"


def _smooth(t):
    t = np.clip(t, 0., 1.)
    return t**3*(10+t*(-15+6*t)), 30*t*t*(1-t)**2


def _box_weight(value, low, high, width):
    wl, dl = _smooth((value-low)/width)
    wh, dh = _smooth((high-value)/width)
    return wl*wh, (dl*wh-wl*dh)/width


def _hermite(t, spacing):
    return (np.array([2*t**3-3*t*t+1, -2*t**3+3*t*t,
                      spacing*(t**3-2*t*t+t), spacing*(t**3-t*t)]),
            np.array([6*t*t-6*t, -6*t*t+6*t,
                      spacing*(3*t*t-4*t+1), spacing*(3*t*t-2*t)])/spacing)


class FourierAuPotential(PeriodicAuPotential):
    def __init__(self, elements, reference_positions=None, *, accelerator="auto",
                 fourier_cache_dir=None, angular_step=.2, chart_bound=2.,
                 chart_blend=.4, height_min=2.8, height_max=8., height_blend=.3,
                 height_step=.05, fourier_nx=8, fourier_ny=18,
                 max_harmonic_x=2, max_harmonic_y=4, **kwargs):
        started = time.perf_counter()
        if accelerator not in ("auto", "numpy", "numba"):
            raise ValueError("Fourier accelerator must be auto, numpy or numba.")
        super().__init__(elements, reference_positions, accelerator="numpy", **kwargs)
        self.accelerator_requested = accelerator
        self.angular_step, self.chart_bound = float(angular_step), float(chart_bound)
        self.chart_blend = float(chart_blend)
        self.height_min, self.height_max = float(height_min), float(height_max)
        self.height_blend, self.height_step = float(height_blend), float(height_step)
        self.fourier_nx, self.fourier_ny = int(fourier_nx), int(fourier_ny)
        if (not np.isfinite([angular_step, chart_bound, chart_blend, height_min,
                            height_max, height_blend, height_step]).all()
            or angular_step <= 0 or chart_bound <= chart_blend or chart_blend <= 0
            or height_max-height_min <= 2*height_blend or height_step <= 0
            or height_min <= 0 or height_blend <= 0):
            raise ValueError("Invalid bounded Fourier interpolation support.")
        angular_intervals = int(round(2*chart_bound/angular_step))
        if not np.isclose(angular_intervals*angular_step, 2*chart_bound):
            raise ValueError("The chart width must be an integer multiple of angular_step.")
        if (self.fourier_nx < 4 or self.fourier_ny < 6
            or self.fourier_nx != fourier_nx or self.fourier_ny != fourier_ny
            or max_harmonic_x < 0 or max_harmonic_y < 0
            or int(max_harmonic_x) != max_harmonic_x or int(max_harmonic_y) != max_harmonic_y
            or self.fourier_nx % 2 or self.fourier_ny % 6
            or self.fourier_nx <= 2*max_harmonic_x or self.fourier_ny <= 2*max_harmonic_y):
            raise ValueError("Fourier sampling must resolve retained harmonics and ABC shifts.")
        self.angular_nodes = np.linspace(-chart_bound, chart_bound, angular_intervals+1)
        self.fourier_z_start = height_min-3*height_step
        self.fourier_z_count = int(np.ceil((height_max-height_min)/height_step))+7
        modes = [(m, n) for m in range(int(max_harmonic_x)+1)
                 for n in range(-int(max_harmonic_y), int(max_harmonic_y)+1)
                 if (m+n) % 2 == 0 and (m > 0 or n >= 0)]
        self.modes = np.asarray(modes, dtype=np.int64)
        self.wavevectors = self.modes*(2*np.pi/self.period)
        self.mode_multiplicity = np.where(np.any(self.modes != 0, axis=1), 2., 1.)
        self.fourier_cache_dir = Path(fourier_cache_dir) if fourier_cache_dir else ROOT/"outputs/fourier-au-check/cache"
        self.fourier_cache_loaded = False
        self.fourier_cache_build_seconds = 0.
        self.last_evaluation_counts = {"full_proxy":0, "blended_proxy":0, "exact_direction_fallback":0}
        self._prepare_fourier_cache()
        self._compiled_fourier_kernel = None
        if accelerator != "numpy" and importlib.util.find_spec("numba"):
            start_import = time.perf_counter()
            from .fourier_au_kernel import energy_and_forces, NUMBA_VERSION
            self.numba_import_seconds = time.perf_counter()-start_import
            self.numba_version = NUMBA_VERSION
            self._compiled_fourier_kernel = energy_and_forces
            if self.hetero_surface is not None:
                self.hetero_surface.accelerator = "numba"
            rows = []
            for element in ("C", "H"):
                pair = []
                for key in ((element, "Au"), ("Au", element)):
                    p = self.parameters[key]
                    pair.append([p.beta, p.decay, p.delta**2, p.epsilon, p.c,
                                 p.d, p.damping_radius, p.c6])
                rows.append(pair)
            self._kernel_parameters = np.asarray(rows, dtype=np.float64)
            self._kernel_grid = self._grid if self._grid is not None else np.zeros((2, 1, 4, 6))
            self._kernel_z_start = self._z_start if self._grid is not None else 0.
            if reference_positions is None:
                reference_positions = np.column_stack([np.arange(len(self.elements))*3.,
                    np.zeros(len(self.elements)), np.full(len(self.elements), 3.4)])
            warmup = time.perf_counter()
            self._evaluate_fourier_compiled(self._positions(reference_positions))
            self.numba_compile_warmup_seconds = time.perf_counter()-warmup
            self.accelerator = "numba"
        elif accelerator == "numba":
            raise ImportError("Fourier acceleration requires the project's Numba dependency.")
        self.preparation_seconds = time.perf_counter()-started

    def _coefficients_at_height(self, z):
        """Finite-cutoff lattice convolution; no Au pair approximation here."""
        nx, ny = self.fourier_nx, self.fourier_ny
        dx, dy = self.period/[nx, ny]
        qx, qy = np.meshgrid(np.arange(-int(np.ceil(16/dx)), int(np.ceil(16/dx))+1),
                            np.arange(-int(np.ceil(16/dy)), int(np.ceil(16/dy))+1), indexing="ij")
        px, py = (qx*dx).ravel(), (qy*dy).ravel()
        lateral = px*px+py*py
        bins = ((qx % nx)*ny+(qy % ny)).ravel()
        shape = (2, len(self.angular_nodes), len(self.angular_nodes), 4, len(self.modes))
        result = np.zeros(shape, dtype=np.complex128)
        aa, bb = np.meshgrid(self.angular_nodes, self.angular_nodes, indexing="ij")
        a, b = aa.ravel(), bb.ravel()
        q = np.sqrt(1+a*a+b*b)
        for layer, height in enumerate(self.layer_heights):
            dz = z-height
            keep = lateral+dz*dz < 256.
            if not np.any(keep):
                continue
            xx, yy, ll, ii = px[keep], py[keep], lateral[keep], bins[keep]
            radius = np.sqrt(ll+dz*dz)
            taper, _ = _taper(radius)
            shift = (layer % 3)*np.array([self.spacing/2, self.period[1]/6])
            phase = np.exp(-1j*(self.wavevectors@shift))
            # Both sites of the rectangular cell; odd m+n modes are absent.
            phase *= 2.
            for code, element in enumerate(("C", "H")):
                p = self.parameters[element, "Au"]
                radial = taper*p.c*np.exp(-p.decay*(radius-p.beta))
                fixed = radial*np.exp(-ll/p.delta**2) if code else None
                for index in range(len(a)):
                    s = (xx*a[index]+yy*b[index]+dz)/q[index]
                    sa = xx/q[index]-s*a[index]/q[index]**2
                    sb = yy/q[index]-s*b[index]/q[index]**2
                    sab = -(xx*b[index]+yy*a[index])/q[index]**3+3*s*a[index]*b[index]/q[index]**4
                    value = radial*np.exp(-(ll+dz*dz-s*s)/p.delta**2)
                    ga, gb = 2*s*sa/p.delta**2, 2*s*sb/p.delta**2
                    derivatives = (value-(fixed if fixed is not None else 0), value*ga, value*gb,
                        value*(ga*gb+2*(sa*sb+s*sab)/p.delta**2))
                    ia, ib = divmod(index, len(self.angular_nodes))
                    for component, values in enumerate(derivatives):
                        field = np.bincount(ii, weights=values, minlength=nx*ny).reshape(nx, ny)
                        transformed = np.fft.fft2(field)/(nx*ny)
                        result[code, ia, ib, component] += transformed[self.modes[:,0] % nx,
                            self.modes[:,1] % ny]*phase
        return result

    def _prepare_fourier_cache(self):
        settings = {"version":ALGORITHM_VERSION, "parameters":self.parameter_sha256,
            "period":self.period.tolist(), "heights":self.layer_heights.tolist(),
            "angular_nodes":self.angular_nodes.tolist(), "z_start":self.fourier_z_start,
            "z_count":self.fourier_z_count, "z_step":self.height_step,
            "xy_sampling":[self.fourier_nx,self.fourier_ny], "modes":self.modes.tolist(),
            "domain":[self.height_min,self.height_max,self.chart_bound,self.chart_blend,self.height_blend]}
        serialized = json.dumps(settings, sort_keys=True)
        identity = hashlib.sha256(serialized.encode()).hexdigest()
        self.fourier_cache_path = self.fourier_cache_dir/f"fourier-{identity}.npz"
        expected = (2,self.fourier_z_count,len(self.angular_nodes),len(self.angular_nodes),4,len(self.modes))
        if self.fourier_cache_path.is_file():
            try:
                with np.load(self.fourier_cache_path, allow_pickle=False) as saved:
                    coeff = saved["coefficients"]
                    if (str(saved["settings"].item()) == serialized and coeff.shape == expected
                        and coeff.dtype == np.dtype(np.complex128) and np.isfinite(coeff).all()):
                        self.coefficients = coeff
                        self.fourier_cache_loaded = True
                        return
            except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile):
                pass
        start = time.perf_counter()
        self.coefficients = np.empty(expected, dtype=np.complex128)
        for iz in range(self.fourier_z_count):
            self.coefficients[:,iz] = self._coefficients_at_height(self.fourier_z_start+iz*self.height_step)
        self.fourier_cache_build_seconds = time.perf_counter()-start
        self.fourier_cache_dir.mkdir(parents=True, exist_ok=True)
        temporary = self.fourier_cache_path.with_name(self.fourier_cache_path.name+f".{os.getpid()}.{time.time_ns()}.tmp")
        try:
            with temporary.open("wb") as output:
                np.savez_compressed(output, coefficients=self.coefficients, settings=serialized)
            os.replace(temporary,self.fourier_cache_path)
        finally:
            if temporary.exists():
                temporary.unlink()

    def _proxy_weight(self, z, normal):
        if abs(normal[2]) < 1e-12:
            return 0., 0., np.zeros(3), 0., 0.
        a,b = normal[:2]/normal[2]
        wz,dz = _box_weight(z,self.height_min,self.height_max,self.height_blend)
        wa,da = _box_weight(a,-self.chart_bound,self.chart_bound,self.chart_blend)
        wb,db = _box_weight(b,-self.chart_bound,self.chart_bound,self.chart_blend)
        grad = np.array([wz*da*wb/normal[2],wz*wa*db/normal[2],
            -wz*(da*wb*a+wa*db*b)/normal[2]])
        # The proxy describes unit normals. A C2 potential weight (including
        # its normal derivative) returns nearly degenerate regularized planes
        # to the complete SAIP sum without a hard-switch force discontinuity.
        wn,dn = _smooth((np.dot(normal,normal)-.999)/.0009)
        weight = wz*wa*wb
        return (float(weight*wn),float(dz*wa*wb*wn),
                grad*wn+weight*dn/.0009*2*normal,float(a),float(b))

    def _proxy(self, code, position, normal, a,b):
        az,bz = (a+self.chart_bound)/self.angular_step,(b+self.chart_bound)/self.angular_step
        iz = (position[2]-self.fourier_z_start)/self.height_step
        ia,ib,kz = int(np.floor(az)),int(np.floor(bz)),int(np.floor(iz))
        ha,da = _hermite(az-ia,self.angular_step)
        hb,db = _hermite(bz-ib,self.angular_step)
        wz,dz = _cubic_weights(np.array([iz-kz]))
        value = np.zeros(len(self.modes),dtype=np.complex128)
        ga,gb,gz = value.copy(),value.copy(),value.copy()
        for j in range(4):
            for k in range(4):
                component = (j//2)+2*(k//2)
                # storage order [value, da, db, da_db]
                for l in range(4):
                    c = self.coefficients[code,kz+l-1,ia+j%2,ib+k%2,component]
                    value += c*ha[j]*hb[k]*wz[0,l]
                    ga += c*da[j]*hb[k]*wz[0,l]
                    gb += c*ha[j]*db[k]*wz[0,l]
                    gz += c*ha[j]*hb[k]*dz[0,l]/self.height_step
        phase = np.exp(1j*(self.wavevectors@(position[:2]%self.period)))*self.mode_multiplicity
        energy = np.sum((phase*value).real)
        gradient = np.array([np.sum((1j*self.wavevectors[:,0]*phase*value).real),
                             np.sum((1j*self.wavevectors[:,1]*phase*value).real),
                             np.sum((phase*gz).real)])
        d_a,d_b = np.sum((phase*ga).real),np.sum((phase*gb).real)
        direction = np.array([d_a/normal[2],d_b/normal[2],-(a*d_a+b*d_b)/normal[2]])
        return float(energy),-gradient,direction

    def _evaluate_fourier_compiled(self,positions,center_weights=None,center_dz=None):
        energy,force,counts = self._compiled_fourier_kernel(positions,self._element_codes,
            self.normal_cutoffs,self._tile,self._kernel_parameters,self._kernel_grid,
            float(self.period[0]),float(self.period[1]),self.grid_dz,self._kernel_z_start,
            self._grid is not None,self.coefficients,self.wavevectors,self.mode_multiplicity,
            self.angular_step,self.chart_bound,self.chart_blend,self.height_step,
            self.fourier_z_start,self.height_min,self.height_max,self.height_blend,
            self._normal_neighbours,self._normal_counts,
            NORMAL_AREA_EPSILON if self.normal_policy == BONDED_NORMAL_POLICY else 0.,center_weights,center_dz)
        self.last_evaluation_counts = dict(zip(("full_proxy","blended_proxy","exact_direction_fallback"),map(int,counts)))
        return self._add_hetero((float(energy),force), positions,center_weights,center_dz)

    def energy_and_forces_eV(self,positions,center_weights=None,center_dz=None):
        positions = self._positions(positions)
        if self._compiled_fourier_kernel is not None:
            return self._evaluate_fourier_compiled(positions,center_weights,center_dz)
        geometry = self.normal_geometry(positions)
        normals = geometry[0]
        n = len(positions)
        radial_weight,radial_derivative = _cache_weight(positions[:,2]) if self._grid is not None else (np.zeros(n),np.zeros(n))
        cached_e,cached_f = np.zeros(n),np.zeros_like(positions)
        active = np.flatnonzero((radial_weight>0) & (self._element_codes >= 0))
        if len(active):
            cached_e[active],cached_f[active] = self._interpolate(positions,active)
        base_e,base_f = np.zeros(n),np.zeros_like(positions)
        force,normal_gradient = np.zeros_like(positions),np.zeros_like(positions)
        total = 0.
        counts = [0,0,0]
        for i in range(n):
            code = self._element_codes[i]
            if code < 0:
                continue
            direction_dependent = code == 0 or geometry[2][i] >= 2
            w,dw,gn,a,b = self._proxy_weight(positions[i,2],normals[i]) if direction_dependent else (1.,0.,np.zeros(3),0.,0.)
            if direction_dependent:
                counts[0 if w == 1 else 1 if w>0 else 2] += 1
            pe,pf,pg = self._proxy(code,positions[i],normals[i],a,b) if direction_dependent and w>0 else (0.,np.zeros(3),np.zeros(3))
            ee,ef,eg = 0.,np.zeros(3),np.zeros(3)
            if radial_weight[i]<1 or (direction_dependent and w<1):
                ids,delta,radius = self._pairs(positions,np.array([i]))
                p,gold = self.parameters[self.elements[i],"Au"],self.parameters["Au",self.elements[i]]
                fixed_e,fixed_f = np.zeros(len(ids)),np.zeros_like(delta)
                if code == 1:
                    zz = np.zeros_like(delta); zz[:,2]=1
                    fixed_e,fixed_f,_ = _anisotropic(delta,radius,zz,p)
                if radial_weight[i]<1:
                    re,rf = _radial(delta,radius,p,gold)
                    base_e[i],base_f[i] = np.sum(re+fixed_e),np.sum(rf+fixed_f,axis=0)
                if direction_dependent and w<1:
                    ae,af,ag = _anisotropic(delta,radius,np.repeat(normals[i,None],len(ids),axis=0),p)
                    ee,ef,eg = np.sum(ae-fixed_e),np.sum(af-fixed_f,axis=0),np.sum(ag,axis=0)
            atomic_energy = (1-radial_weight[i])*base_e[i]+radial_weight[i]*cached_e[i]+w*pe+(1-w)*ee
            total += atomic_energy if center_weights is None else center_weights[i]*atomic_energy
            force[i] += (1-radial_weight[i])*base_f[i]+radial_weight[i]*cached_f[i]+w*pf+(1-w)*ef
            force[i,2] -= radial_derivative[i]*(cached_e[i]-base_e[i])+dw*(pe-ee)
            normal_gradient[i] = w*pg+(1-w)*eg+gn*(pe-ee)
            if center_weights is not None:
                force[i] *= center_weights[i]
                normal_gradient[i] *= center_weights[i]
                force[i,2] -= center_dz[i]*atomic_energy
        _add_normal_chain_forces(force,normal_gradient,geometry)
        self.last_evaluation_counts = dict(zip(("full_proxy","blended_proxy","exact_direction_fallback"),counts))
        return self._add_hetero((float(total),force), positions,center_weights,center_dz)

    def metadata(self):
        result = super().metadata()
        result.update({"model":"bounded Fourier proxy of frozen six-layer Au(111), SAIP CHAu",
            "algorithm_version":ALGORITHM_VERSION,"evaluation":"radial/H cache + bounded full-direction Fourier/Hermite proxy with exact conservative fallback",
            "fourier_modes":self.modes.tolist(),"fourier_mode_count":len(self.modes),
            "direction_representation":"a=nx/nz, b=ny/nz; bicubic Hermite of full SAIP normal dependence",
            "proxy_height_support_angstrom":[self.height_min,self.height_max],
            "proxy_height_blend_angstrom":self.height_blend,"proxy_chart_bound":self.chart_bound,
            "proxy_chart_blend":self.chart_blend,"proxy_angular_step":self.angular_step,
            "proxy_height_step_angstrom":self.height_step,"fourier_xy_sampling":[self.fourier_nx,self.fourier_ny],
            "fourier_cache_shape":list(self.coefficients.shape),"fourier_cache_file":str(self.fourier_cache_path),
            "fourier_cache_loaded":self.fourier_cache_loaded,"fourier_cache_build_seconds":self.fourier_cache_build_seconds,
            "last_evaluation_counts":self.last_evaluation_counts.copy(),
            "fallback":"per-atom complete 16 A direction sum outside bounded proxy; smooth potential blending including position and normal weight derivatives"})
        result["proxy_regularized_normal_support"] = {
            "normal_squared_magnitude_exact_below":.999,
            "normal_squared_magnitude_proxy_above":.9999,
            "blend":"C2 scalar potential weight with full normal derivative"}
        if self.hetero_surface is not None:
            result["model"] = "bounded Fourier CHAu proxy + provisional X-Au game LJ"
        return result


def create_fourier_surface(atoms,reference_positions=None,**kwargs):
    return FourierAuPotential(atoms,reference_positions,**kwargs)
