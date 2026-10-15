"""Three-chart Fourier v2 for the native game and independent reference.

The normalized C2 partition covers every nondegenerate unit normal. All chart
and height derivatives act on the scalar energy before normal-chain forces.
The native game uses a compiled, fused-phase kernel; accelerator="numpy"
retains the independent Python reference for energy and force comparisons.
"""
import hashlib
import importlib.util
import json
import time
import numpy as np
from .fourier_au import FourierAuPotential, _box_weight, _smooth
from .periodic_au import _taper
from .periodic_au import NORMAL_AREA_EPSILON, BONDED_NORMAL_POLICY

AXES = ((0, 1, 2), (1, 2, 0), (2, 0, 1))


class FourierAtlasPotential(FourierAuPotential):
    def __init__(self, *args, accelerator="auto", **kwargs):
        if accelerator not in ("auto", "numpy", "numba"):
            raise ValueError("Atlas accelerator must be auto, numpy or numba.")
        # Build the unchanged v2 tables without installing the v1 kernel.
        kwargs.update(accelerator="numpy")
        kwargs.setdefault("angular_step", .2)
        kwargs.setdefault("chart_bound", 1.2)
        kwargs.setdefault("chart_blend", .2)
        kwargs.setdefault("height_min", 2.4)
        kwargs.setdefault("height_max", 11.2)
        kwargs.setdefault("height_step", .1)
        super().__init__(*args, **kwargs)
        self.accelerator_requested = accelerator
        self._phase_mean = np.asarray([(m % 2 == 0 and n % 3 == 0)
                                      for m, n in self.modes], dtype=float)
        self._unit_weights = np.ones(len(self.elements))
        self._zero_dz = np.zeros(len(self.elements))
        if accelerator != "numpy" and importlib.util.find_spec("numba"):
            started = time.perf_counter()
            from .fourier_atlas_kernel import energy_and_forces, NUMBA_VERSION
            self.numba_import_seconds = time.perf_counter()-started
            self.numba_version = NUMBA_VERSION
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
            self._compiled_fourier_kernel = energy_and_forces
            if self.hetero_surface is not None:
                self.hetero_surface.accelerator = "numba"
            warmup = np.column_stack([np.arange(len(self.elements))*3.,
                                      np.zeros(len(self.elements)), np.full(len(self.elements), 3.4)])
            started = time.perf_counter()
            self._evaluate_fourier_compiled(warmup, scale=1.5)
            self.numba_compile_warmup_seconds = time.perf_counter()-started
            self.accelerator = "numba"
            self.preparation_seconds += self.numba_import_seconds+self.numba_compile_warmup_seconds
        elif accelerator == "numba":
            raise ImportError("Atlas acceleration requires the project's Numba dependency.")

    def _evaluate_fourier_compiled(self, positions, center_weights=None, center_dz=None, *, scale=1.):
        weights = self._unit_weights if center_weights is None else center_weights
        dz = self._zero_dz if center_weights is None else center_dz
        energy, force, counts = self._compiled_fourier_kernel(
            positions, self._element_codes, self.normal_cutoffs, self._tile,
            self._kernel_parameters, self._kernel_grid, float(self.period[0]), float(self.period[1]),
            self.grid_dz, self._kernel_z_start, self._grid is not None,
            self.atlas_coefficients, self.wavevectors, self.mode_multiplicity,
            self.angular_step, self.chart_bound, self.chart_blend, self.height_step,
            self.fourier_z_start, self.height_min, self.height_max, self.height_blend,
            self._normal_neighbours, self._normal_counts,
            NORMAL_AREA_EPSILON if self.normal_policy == BONDED_NORMAL_POLICY else 0.,
            self._phase_mean, float(scale), weights, dz)
        self.last_evaluation_counts = dict(zip(
            ("full_proxy", "blended_proxy", "exact_direction_fallback"), map(int, counts)))
        if self.hetero_surface is None or scale == 1.:
            return self._add_hetero((float(energy), force), positions, center_weights, center_dz)
        # X-Au remains the existing compiled radial model. Apply the same
        # atomic height gate and exact six translations as the reference.
        t = np.clip((positions[:, 2]-2.)/.6, 0., 1.)
        factor = 1.+(scale-1.)*t**3*(10+t*(-15+6*t))
        derivative = (scale-1.)*30*t*t*(t-1)**2/.6
        for i in range(2):
            for j in range(3):
                w = (1+5*factor)/6 if i == j == 0 else (1-factor)/6
                d = 5*derivative/6 if i == j == 0 else -derivative/6
                energy, force = self._add_hetero((energy, force),
                    positions+[i*self.period[0]/2, j*self.period[1]/3, 0.],
                    weights*w, dz*w+weights*d)
        return float(energy), force

    def energy_and_forces_corrugated_eV(self, positions, scale):
        return self._evaluate_fourier_compiled(self._positions(positions), scale=scale)

    def _coefficients_at_height(self, z):
        nx, ny = self.fourier_nx, self.fourier_ny
        dx, dy = self.period / [nx, ny]
        qx, qy = np.meshgrid(np.arange(-int(np.ceil(16/dx)), int(np.ceil(16/dx))+1),
                            np.arange(-int(np.ceil(16/dy)), int(np.ceil(16/dy))+1), indexing="ij")
        px, py = (qx*dx).ravel(), (qy*dy).ravel()
        lateral = px*px+py*py
        bins = ((qx % nx)*ny+(qy % ny)).ravel()
        result = np.zeros((2, 3, len(self.angular_nodes), len(self.angular_nodes), 4, len(self.modes)), complex)
        aa, bb = np.meshgrid(self.angular_nodes, self.angular_nodes, indexing="ij")
        a, b = aa.ravel(), bb.ravel()
        q = np.sqrt(1+a*a+b*b)
        for layer, height in enumerate(self.layer_heights):
            dz = z-height
            keep = lateral+dz*dz < 256
            if not np.any(keep):
                continue
            xx, yy, ll, ii = px[keep], py[keep], lateral[keep], bins[keep]
            rr = np.sqrt(ll+dz*dz)
            taper, _ = _taper(rr)
            shift = (layer % 3)*np.array([self.spacing/2, self.period[1]/6])
            phase = 2*np.exp(-1j*(self.wavevectors@shift))
            components = (xx, yy, np.full_like(xx, dz))
            for code, element in enumerate(("C", "H")):
                p = self.parameters[element, "Au"]
                radial = taper*p.c*np.exp(-p.decay*(rr-p.beta))
                fixed = radial*np.exp(-ll/p.delta**2) if code else 0
                for chart, axes in enumerate(AXES):
                    ra, rb, rc = (components[k] for k in axes)
                    for index in range(len(a)):
                        s = (ra*a[index]+rb*b[index]+rc)/q[index]
                        sa = ra/q[index]-s*a[index]/q[index]**2
                        sb = rb/q[index]-s*b[index]/q[index]**2
                        sab = -(ra*b[index]+rb*a[index])/q[index]**3+3*s*a[index]*b[index]/q[index]**4
                        value = radial*np.exp(-(ll+dz*dz-s*s)/p.delta**2)
                        ga, gb = 2*s*sa/p.delta**2, 2*s*sb/p.delta**2
                        derivatives = (value-fixed, value*ga, value*gb,
                            value*(ga*gb+2*(sa*sb+s*sab)/p.delta**2))
                        ia, ib = divmod(index, len(self.angular_nodes))
                        for component, values in enumerate(derivatives):
                            field = np.bincount(ii, weights=values, minlength=nx*ny).reshape(nx, ny)
                            transformed = np.fft.fft2(field)/(nx*ny)
                            result[code, chart, ia, ib, component] += transformed[self.modes[:,0] % nx,
                                self.modes[:,1] % ny]*phase
        return result

    def _prepare_fourier_cache(self):
        settings = dict(version="three-chart-hermite-v2", parameters=self.parameter_sha256,
            period=self.period.tolist(), layers=self.layer_heights.tolist(),
            angular_nodes=self.angular_nodes.tolist(), zstart=self.fourier_z_start,
            nz=self.fourier_z_count, dz=self.height_step, modes=self.modes.tolist(),
            sampling=[self.fourier_nx, self.fourier_ny])
        serialized = json.dumps(settings, sort_keys=True)
        key = hashlib.sha256(serialized.encode()).hexdigest()
        self.fourier_cache_dir.mkdir(parents=True, exist_ok=True)
        self.fourier_cache_path = self.fourier_cache_dir/f"atlas-{key}.npz"
        shape = (2, 3, self.fourier_z_count, len(self.angular_nodes), len(self.angular_nodes), 4, len(self.modes))
        if self.fourier_cache_path.exists():
            with np.load(self.fourier_cache_path, allow_pickle=False) as saved:
                assert str(saved['settings'].item()) == serialized
                self.atlas_coefficients = saved['coefficients']
                assert self.atlas_coefficients.shape == shape
                assert np.isfinite(self.atlas_coefficients).all()
            self.fourier_cache_loaded = True
        else:
            start = time.perf_counter()
            self.atlas_coefficients = np.empty(shape, complex)
            for iz in range(self.fourier_z_count):
                self.atlas_coefficients[:, :, iz] = self._coefficients_at_height(self.fourier_z_start+iz*self.height_step)
                if iz % 10 == 0:
                    print(f"atlas height {iz+1}/{self.fourier_z_count}", flush=True)
            np.savez_compressed(self.fourier_cache_path, coefficients=self.atlas_coefficients, settings=serialized)
            self.fourier_cache_build_seconds = time.perf_counter()-start
        self.coefficients = self.atlas_coefficients[:, 0]

    def _proxy_weight(self, z, normal):
        wz, dz = _box_weight(z, self.height_min, self.height_max, self.height_blend)
        wn, dn = _smooth((np.dot(normal, normal)-.999)/.0009)
        return float(wz*wn), float(dz*wn), wz*dn/.0009*2*normal, 0., 0.

    def _proxy(self, code, position, normal, a, b):
        charts = []
        for chart, axes in enumerate(AXES):
            n = normal[list(axes)]
            if abs(n[2]) < 1e-12:
                continue
            a, b = n[:2]/n[2]
            wa, da = _box_weight(a, -self.chart_bound, self.chart_bound, self.chart_blend)
            wb, db = _box_weight(b, -self.chart_bound, self.chart_bound, self.chart_blend)
            w = wa*wb
            if w == 0:
                continue
            gradient = np.array([da*wb, wa*db, -(a*da*wb+b*wa*db)])/n[2]
            self.coefficients = self.atlas_coefficients[:, chart]
            e, f, ng = super()._proxy(code, position, n, a, b)
            global_g, global_ng = np.empty(3), np.empty(3)
            global_g[list(axes)], global_ng[list(axes)] = gradient, ng
            charts.append((w, global_g, e, f, global_ng))
        total = sum(c[0] for c in charts)
        assert total > 0
        e = sum(w*ce for w, _, ce, _, _ in charts)/total
        f = sum(w*cf for w, _, _, cf, _ in charts)/total
        ng = sum(w*cg+gw*(ce-e) for w, gw, ce, _, cg in charts)/total
        return e, f, ng

    def metadata(self):
        result = super().metadata()
        result.update(algorithm_version="three-chart-hermite-v2",
            direction_representation="three projective normal charts; normalized C2 scalar potential partition",
            chart_axes=[list(a) for a in AXES],
            fourier_cache_shape=list(self.atlas_coefficients.shape),
            evaluation=("compiled v2 atlas with fused six-phase Fourier modes" if self.accelerator == "numba"
                        else "Python reference v2 atlas")+"; exact conservative height/degeneracy fallback")
        return result
