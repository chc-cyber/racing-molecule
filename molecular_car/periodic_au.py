"""Frozen, periodic six-layer Au(111) environment derived from CHAu.ILP.

The molecule retains its complete unwrapped coordinates. Only the substrate
query coordinates are reduced modulo its rectangular lattice repeat. Every
physical Au site inside the original 16 A cutoff is summed once; the molecular
normal neighbours are never made periodic. The default distance policy follows
frozen SAIP; the explicitly selected bonded policy is a conservative game
adaptation for folded, sp and sp3 molecules. Neither models electrons or a
thermally moving gold substrate.

The public game API uses kcal/mol and kcal/mol/A; the explicit eV API supports
comparison with LAMMPS ``units metal``. C/H normals and their coordinate chain
forces follow the LAMMPS ILP/SAIP definition in the reference policy, including
fixed-z normals with at most one short-range neighbour. The bonded policy
retains that directional potential with immutable normal classes and a
regularized local plane; it does not assert SAIP transferability to all C/H.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

import numpy as np
from .element_models import SUPPORTED_RACE_ELEMENTS

ROOT = Path(__file__).resolve().parents[1]
EV_TO_KCAL_MOL = 23.0605478306
CUTOFF = 16.0
ALGORITHM_VERSION = "periodic-saip-2"
CACHE_Z_MIN = 2.0
CACHE_Z_MAX = 12.0
CACHE_BLEND_WIDTH = .3
REFERENCE_NORMAL_POLICY = "distance-saip-reference-v1"
BONDED_NORMAL_POLICY = "bonded-saip-game-v1"
NORMAL_AREA_EPSILON = 1e-6  # A^2; only the explicitly selected game adaptation.


@dataclass(frozen=True)
class Parameters:
    beta: float
    alpha: float
    delta: float
    epsilon: float
    c: float
    d: float
    sr: float
    reff: float
    c6: float
    normal_cutoff: float

    @property
    def decay(self):
        return self.alpha / self.beta

    @property
    def damping_radius(self):
        return self.sr * self.reff


def parameter_file() -> Path:
    installation = os.environ.get("MOLECULAR_CAR_LAMMPS_HOME")
    config = ROOT / ".local/lammps.json"
    if not installation and config.is_file():
        installation = json.loads(config.read_text(encoding="utf-8-sig")).get("installation")
    if not installation:
        raise FileNotFoundError("Configure LAMMPS before preparing a periodic Au surface.")
    return Path(installation) / "Potentials/CHAu.ILP"


def read_parameters(path: Path):
    result = {}
    payload = path.read_bytes()
    for line in payload.decode("utf-8").splitlines():
        row = line.split("#", 1)[0].split()
        if not row:
            continue
        if len(row) != 13:
            raise ValueError("Expected the eleven CHAu.ILP numeric parameters.")
        beta, alpha, delta, epsilon, c, d, sr, reff, c6, scale, normal_cutoff = map(float, row[2:])
        # The file stores meV, with the per-row energy scale S.
        result[tuple(row[:2])] = Parameters(beta, alpha, delta, epsilon*scale*.001,
            c*scale*.001, d, sr, reff, c6*scale*.001, normal_cutoff)
    for a in ("C", "H"):
        for b in ("C", "H", "Au"):
            if (a, b) not in result:
                raise ValueError(f"Missing SAIP coefficient {a}/{b}.")
        if ("Au", a) not in result:
            raise ValueError(f"Missing SAIP coefficient Au/{a}.")
        p, q = result[a, "Au"], result["Au", a]
        if not np.allclose([p.c6, p.d, p.damping_radius], [q.c6, q.d, q.damping_radius], rtol=0, atol=1e-12):
            raise ValueError("The periodic model requires symmetric Au/CH dispersion parameters.")
    return result, hashlib.sha256(payload).hexdigest()


def reference_layer_heights(lattice_constant=4.08):
    if float(lattice_constant) != 4.08:
        return np.arange(-5, 1)*lattice_constant/np.sqrt(3), "ideal FCC(111), caller-selected lattice constant"
    source = Path(__file__).resolve().parent / "data/au111_reference.json"
    data = json.loads(source.read_text(encoding="utf-8"))
    heights = np.asarray(data["layer_heights_angstrom"], dtype=float)
    if (data.get("schema") != "frozen-au111-reference-v1"
        or data.get("lattice_constant_angstrom") != 4.08 or heights.shape != (6,)
        or not np.isfinite(heights).all() or not np.all(np.diff(heights) > 0)
        or heights[-1] != 0.):
        raise ValueError("Invalid bundled frozen Au(111) reference.")
    return heights, str(source)


def _bonded_normal_topology(elements, bonds, reference_positions):
    """Immutable graph and normal classes for the game SAIP adaptation.

    This is not the original distance-selected SAIP definition. No nearest-three
    selection is made: nonplanar four-coordinate sites and linear sp sites use
    the surface axis, while valid two/three-neighbour planes retain all bonds.
    The initial collinearity classification is never revisited during dynamics.
    """
    if bonds is None or reference_positions is None:
        raise ValueError("Bonded game normals require explicit bonds and reference positions.")
    adjacency = [[] for _ in elements]
    edges = []
    seen = set()
    for bond in bonds:
        a, b, order = bond.get("a"), bond.get("b"), bond.get("order", 1.)
        if (isinstance(a, bool) or isinstance(b, bool) or not isinstance(a, (int, np.integer))
            or not isinstance(b, (int, np.integer)) or not 1 <= a <= len(elements)
            or not 1 <= b <= len(elements) or a == b):
            raise ValueError("Normal topology requires distinct one-based atom indices.")
        order = float(order)
        if not np.isfinite(order) or order not in (1., 1.5, 2., 3.):
            raise ValueError("Normal topology requires single, aromatic, double or triple bonds.")
        edge = tuple(sorted((int(a)-1, int(b)-1)))
        if edge in seen:
            raise ValueError("Duplicate bond in normal topology.")
        seen.add(edge)
        edges.append((*edge, order))
        adjacency[edge[0]].append((edge[1], order))
        adjacency[edge[1]].append((edge[0], order))
    neighbours = np.full((len(elements), 3), -1, dtype=np.int64)
    counts = np.zeros(len(elements), dtype=np.int64)
    modes = []
    for i, entries in enumerate(adjacency):
        entries.sort()
        degree = len(entries)
        if elements[i] not in ("C", "H"):
            # This center has an isotropic X-Au interaction. Its coordinates
            # still contribute to neighbouring C/H planes and their chain forces.
            modes.append("isotropic-hetero-au-game")
            continue
        if (elements[i] == "H" and (degree > 1 or any(order != 1. for _, order in entries))) or degree > 4:
            raise ValueError("Normal topology is incompatible with a C/H covalent molecule.")
        if degree <= 1:
            modes.append("surface-axis-terminal")
            continue
        if degree == 4:
            modes.append("surface-axis-four-coordinate-game")
            continue
        if any(order == 3. for _, order in entries) or (degree == 2 and sum(o for _, o in entries) >= 4.):
            modes.append("surface-axis-linear-sp-game")
            continue
        ids = [j for j, _ in entries]
        pivot = ids[2] if degree == 3 else i
        first = reference_positions[ids[0]]-reference_positions[pivot]
        second = reference_positions[ids[1]]-reference_positions[pivot]
        scale = np.linalg.norm(first)*np.linalg.norm(second)
        if scale == 0 or np.linalg.norm(np.cross(first, second)) <= 1e-5*scale:
            modes.append("surface-axis-initial-collinear-game")
            continue
        counts[i] = degree
        neighbours[i, :degree] = ids
        modes.append("bonded-regularized-local-plane-game")
    digest = hashlib.sha256(json.dumps(sorted(edges), separators=(",", ":")).encode()).hexdigest()
    neighbours.setflags(write=False)
    counts.setflags(write=False)
    return neighbours, counts, tuple(modes), digest


def _bonded_normal_geometry(positions, neighbours, counts, epsilon=NORMAL_AREA_EPSILON):
    normals = np.zeros_like(positions)
    normals[:, 2] = 1.
    first, second = np.zeros_like(positions), np.zeros_like(positions)
    length = np.ones(len(positions))
    active = np.flatnonzero(counts >= 2)
    pivot = np.where(counts[active] == 3, neighbours[active, 2], active)
    first[active] = positions[neighbours[active, 0]]-positions[pivot]
    second[active] = positions[neighbours[active, 1]]-positions[pivot]
    raw = np.cross(first[active], second[active])
    length[active] = np.sqrt(np.sum(raw*raw, axis=1)+epsilon*epsilon)
    normals[active] = raw/length[active, None]
    return normals, neighbours, counts, first, second, length


def _normal_geometry(positions, cutoffs):
    displacement = positions[None, :, :] - positions[:, None, :]
    square = np.einsum("ijk,ijk->ij", displacement, displacement)
    mask = (square > 0) & (square < cutoffs**2)
    counts = mask.sum(axis=1)
    if np.any(counts > 3):
        atom = int(np.flatnonzero(counts > 3)[0])
        raise ValueError(f"SAIP atom {atom+1} has more than three normal neighbours; this geometry is unsupported.")
    neighbours = np.full((len(positions), 3), -1, dtype=int)
    rows, columns = np.nonzero(mask)
    offset = np.repeat(np.cumsum(counts)-counts, counts)
    neighbours[rows, np.arange(len(rows))-offset] = columns
    normals = np.zeros_like(positions)
    normals[:, 2] = 1
    first, second = np.zeros_like(positions), np.zeros_like(positions)
    active = np.flatnonzero(counts >= 2)
    pivot = np.where(counts[active] == 3, neighbours[active, 2], active)
    first[active] = positions[neighbours[active, 0]] - positions[pivot]
    second[active] = positions[neighbours[active, 1]] - positions[pivot]
    raw = np.cross(first[active], second[active])
    length = np.ones(len(positions))
    length[active] = np.linalg.norm(raw, axis=1)
    if np.any(length[active] < 1e-12):
        raise ValueError("A SAIP local plane has zero normal magnitude.")
    normals[active] = raw / length[active, None]
    return normals, neighbours, counts, first, second, length


def _add_normal_chain_forces(forces, normal_gradient, geometry):
    normal, neighbours, counts, first, second, length = geometry
    active = np.flatnonzero(counts >= 2)
    gradient = normal_gradient[active]
    adjoint = (gradient - normal[active]*np.sum(gradient*normal[active], axis=1)[:, None])/length[active, None]
    ga = np.cross(second[active], adjoint)
    gb = np.cross(adjoint, first[active])
    np.add.at(forces, neighbours[active, 0], -ga)
    np.add.at(forces, neighbours[active, 1], -gb)
    pivot = np.where(counts[active] == 3, neighbours[active, 2], active)
    np.add.at(forces, pivot, ga+gb)


def _taper(radius):
    x = radius/CUTOFF
    return 1+x**4*(-35+x*(84+x*(-70+20*x))), 140*x**3*(x-1)**3/CUTOFF


def _radial(delta, radius, p, gold):
    taper, derivative = _taper(radius)
    organic_rep = .5*p.epsilon*np.exp(-p.decay*(radius-p.beta))
    gold_rep = (.5*gold.epsilon+gold.c)*np.exp(-gold.decay*(radius-gold.beta))
    # Stable logistic avoids overflowing at very short separations.
    argument = p.d*(radius/p.damping_radius-1)
    damping = np.exp(-np.logaddexp(0, -argument))
    dispersion = -p.c6*radius**-6*damping
    raw = organic_rep + gold_rep + dispersion
    radial_derivative = -p.decay*organic_rep-gold.decay*gold_rep
    radial_derivative += p.c6*radius**-7*damping*(6-radius*p.d/p.damping_radius*(1-damping))
    energy = taper*raw
    force = -(derivative*raw+taper*radial_derivative)[:, None]*delta/radius[:, None]
    return energy, force


def _anisotropic(delta, radius, normal, p):
    taper, derivative = _taper(radius)
    projection = np.sum(delta*normal, axis=1)
    perpendicular = delta-projection[:, None]*normal
    # This projection form also defines a conservative, continuous extension
    # when the game's regularized plane vector has magnitude below one.
    rho_square = np.sum(delta*delta, axis=1)-projection*projection
    exponential = p.c*np.exp(-p.decay*(radius-p.beta)-rho_square/p.delta**2)
    energy = taper*exponential
    force = -(derivative-p.decay*taper)[:, None]*exponential[:, None]*delta/radius[:, None]
    force += 2*energy[:, None]*perpendicular/p.delta**2
    normal_gradient = 2*energy[:, None]*projection[:, None]*delta/p.delta**2
    return energy, force, normal_gradient


def _cubic_weights(t):
    """Catmull-Rom potential interpolation and its own analytic derivative."""
    square, cube = t*t, t*t*t
    weights = np.stack([-.5*t+square-.5*cube, 1-2.5*square+1.5*cube,
                        .5*t+2*square-1.5*cube, -.5*square+.5*cube], axis=1)
    derivative = np.stack([-.5+2*t-1.5*square, -5*t+4.5*square,
                           .5+4*t-4.5*square, -t+1.5*square], axis=1)
    return weights, derivative


def _cache_weight(z):
    # The fallback is another conservative potential. Include the derivative
    # of this weight in the force, rather than blending force fields alone.
    low = np.clip((z-CACHE_Z_MIN)/CACHE_BLEND_WIDTH, 0, 1)
    high = np.clip((CACHE_Z_MAX-z)/CACHE_BLEND_WIDTH, 0, 1)
    smooth = lambda t: t**3*(10+t*(-15+6*t))
    deriv = lambda t: 30*t*t*(t-1)**2/CACHE_BLEND_WIDTH
    return smooth(low)*smooth(high), deriv(low)*smooth(high)-smooth(low)*deriv(high)


class PeriodicAuPotential:
    """Periodic frozen Au surface; public forces have the same units as fields.py."""

    def __init__(self, elements, reference_positions=None, *, use_cache=True,
                 parameter_path=None, layer_heights=None, lattice_constant=4.08,
                 cache_dir=None, grid_nx=32, grid_ny=60, grid_dz=.025,
                 accelerator="auto", bonds=None, normal_policy=REFERENCE_NORMAL_POLICY):
        start = time.perf_counter()
        if accelerator not in ("auto", "numpy", "numba"):
            raise ValueError("Periodic Au accelerator must be 'auto', 'numpy' or 'numba'.")
        self.accelerator_requested = accelerator
        self.accelerator = "numpy"
        self.numba_version = None
        self.numba_import_seconds = 0.
        self.numba_compile_warmup_seconds = 0.
        self._compiled_kernel = None
        self.elements = tuple(a["element"] if isinstance(a, dict) else str(a) for a in elements)
        self.has_hetero = any(e not in ("C", "H") for e in self.elements)
        if not self.elements or (self.has_hetero and normal_policy != BONDED_NORMAL_POLICY):
            raise ValueError("Periodic Au supports C/H molecules only; other elements have no CHAu SAIP parameters.")
        if any(e not in SUPPORTED_RACE_ELEMENTS for e in self.elements):
            raise ValueError("Periodic Au game extension has no parameters for these elements.")
        if self.has_hetero and not all(isinstance(a, dict) for a in elements):
            raise ValueError("Hetero-Au game surface requires validated atom-type UFF vdW parameters.")
        self.parameter_path = Path(parameter_path) if parameter_path else parameter_file()
        self.parameters, self.parameter_sha256 = read_parameters(self.parameter_path)
        self.normal_cutoffs = np.array([[self.parameters[a, b].normal_cutoff if (a, b) in self.parameters else 0.
                                        for b in self.elements] for a in self.elements])
        self._element_codes = np.array([0 if a == "C" else 1 if a == "H" else -1 for a in self.elements])
        if normal_policy not in (REFERENCE_NORMAL_POLICY, BONDED_NORMAL_POLICY):
            raise ValueError("Unknown SAIP normal policy.")
        self.normal_policy = normal_policy
        self._normal_neighbours = np.full((len(self.elements), 3), -1, dtype=np.int64)
        self._normal_counts = np.zeros(len(self.elements), dtype=np.int64)
        self.normal_modes = None
        self.normal_topology_sha256 = None
        if normal_policy == BONDED_NORMAL_POLICY:
            reference = self._positions(reference_positions) if reference_positions is not None else None
            (self._normal_neighbours, self._normal_counts, self.normal_modes,
             self.normal_topology_sha256) = _bonded_normal_topology(self.elements, bonds, reference)
        if reference_positions is not None:
            self.normal_geometry(self._positions(reference_positions))
        self.spacing = float(lattice_constant)/np.sqrt(2)
        self.period = np.array([self.spacing, np.sqrt(3)*self.spacing])
        default_heights, self.reference_source = reference_layer_heights(lattice_constant)
        self.layer_heights = np.asarray(default_heights if layer_heights is None else layer_heights, dtype=float)
        if self.layer_heights.shape != (6,) or not np.isfinite(self.layer_heights).all():
            raise ValueError("The frozen reference must provide six finite layer heights, bottom to top.")
        if layer_heights is not None:
            self.reference_source = "explicit caller-provided six layer heights"
        self._tile = self._build_tile()
        self._grid = None
        self.use_cache = bool(use_cache)
        self.cache_dir = Path(cache_dir) if cache_dir else ROOT/"outputs/periodic-au-check/cache"
        self.grid_nx, self.grid_ny, self.grid_dz = int(grid_nx), int(grid_ny), float(grid_dz)
        if self.grid_nx < 4 or self.grid_nx % 2 or self.grid_ny < 6 or self.grid_ny % 6:
            raise ValueError("Cache nx must be even and >=4; ny must be divisible by six and >=6.")
        if not np.isfinite(self.grid_dz) or not 0 < self.grid_dz <= .25:
            raise ValueError("Cache vertical spacing must be in (0, .25] angstrom.")
        self.cache_build_seconds = 0.
        self.cache_loaded = False
        self.cache_path = None
        if self.use_cache:
            self._prepare_cache()
        self.hetero_surface = None
        if self.has_hetero:
            from .hetero_au import HeteroAuPotential
            self.hetero_surface = HeteroAuPotential(elements, self)
        self._prepare_accelerator(reference_positions)
        self.preparation_seconds = time.perf_counter()-start

    def _positions(self, positions):
        # A fixed contiguous float64 ABI also prevents a strided view from
        # triggering another JIT specialization during the first MD step.
        values = np.ascontiguousarray(positions, dtype=float)
        if values.shape != (len(self.elements), 3) or not np.isfinite(values).all():
            raise ValueError("Provide finite complete unwrapped molecular coordinates, shape N by 3.")
        return values

    def normal_geometry(self, positions):
        """Use the preparation-time policy on complete unwrapped coordinates."""
        if self.normal_policy == BONDED_NORMAL_POLICY:
            return _bonded_normal_geometry(positions, self._normal_neighbours, self._normal_counts)
        return _normal_geometry(positions, self.normal_cutoffs)

    def _prepare_accelerator(self, reference_positions):
        if self.accelerator_requested == "numpy":
            return
        available = importlib.util.find_spec("numba") is not None
        if not available:
            if self.accelerator_requested == "numba":
                raise ImportError("The numba accelerator requires the project's optional Numba dependency.")
            return
        started = time.perf_counter()
        from .periodic_au_kernel import energy_and_forces, NUMBA_VERSION
        self.numba_import_seconds = time.perf_counter()-started
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
        self._compiled_kernel = energy_and_forces
        if self.hetero_surface is not None:
            self.hetero_surface.accelerator = "numba"
        self.numba_version = NUMBA_VERSION
        if reference_positions is None:
            # Independent isolated atoms provide a valid warmup geometry even
            # when only elements were supplied to construct a generic field.
            reference_positions = np.column_stack([np.arange(len(self.elements))*3.,
                                                   np.zeros(len(self.elements)),
                                                   np.full(len(self.elements), 3.4)])
        started = time.perf_counter()
        self._evaluate_compiled(self._positions(reference_positions))
        self.numba_compile_warmup_seconds = time.perf_counter()-started
        self.accelerator = "numba"

    def _evaluate_compiled(self, positions, center_weights=None, center_dz=None):
        result = self._compiled_kernel(positions, self._element_codes, self.normal_cutoffs,
            self._tile, self._kernel_parameters, self._kernel_grid,
            float(self.period[0]), float(self.period[1]), self.grid_dz,
            self._kernel_z_start, self._grid is not None, self._normal_neighbours,
            self._normal_counts, NORMAL_AREA_EPSILON if self.normal_policy == BONDED_NORMAL_POLICY else 0., center_weights, center_dz)
        return self._add_hetero(result, positions, center_weights, center_dz)

    def _add_hetero(self, result, positions, center_weights=None, center_dz=None):
        if self.hetero_surface is None:
            return result
        energy, forces = result
        extra_energy, extra_forces = self.hetero_surface.energy_and_forces_eV(positions, center_weights, center_dz)
        return float(energy+extra_energy), forces+extra_forces

    def _build_tile(self):
        # Include every lattice site which could be within cutoff of the
        # rectangular query cell; images are distinct physical Au sites.
        nx, ny = np.ceil(CUTOFF/self.period).astype(int)+1
        ix, iy = np.meshgrid(np.arange(-nx, nx+1), np.arange(-ny, ny+1), indexing="ij")
        xy = np.stack([ix.ravel()*self.period[0], iy.ravel()*self.period[1]], axis=1)
        layers = []
        for k, height in enumerate(self.layer_heights):
            shift = (k % 3)*np.array([self.spacing/2, self.period[1]/6])
            for basis in ([0., 0.], [self.spacing/2, self.period[1]/2]):
                sites = xy + np.mod(shift+basis, self.period)
                nearest = np.maximum(np.maximum(-sites, sites-self.period), 0)
                keep = np.sum(nearest**2, axis=1) <= CUTOFF**2
                layers.append(np.column_stack([sites[keep], np.full(np.count_nonzero(keep), height)]))
        return np.concatenate(layers)

    def _pairs(self, positions, active):
        local = positions[active].copy()
        local[:, :2] %= self.period
        delta = local[:, None, :]-self._tile[None, :, :]
        square = np.einsum("ijk,ijk->ij", delta, delta)
        ia, ja = np.nonzero(square < CUTOFF**2)
        radius = np.sqrt(square[ia, ja])
        if np.any(radius < 1e-8):
            raise ValueError("A molecular atom coincides with a gold lattice site.")
        return active[ia], delta[ia, ja], radius

    def _prepare_cache(self):
        # Cache two scalar potentials. Carbon retains its entire direction-
        # dependent term as an exact lattice sum. Fixed-z hydrogen caches the
        # whole potential; hydrogen with a real local normal gets an exact
        # direction-dependent correction, including its normal chain force.
        settings = {"version":ALGORITHM_VERSION,"parameters":self.parameter_sha256,
                    "period":self.period.tolist(),"heights":self.layer_heights.tolist(),
                    "nx":self.grid_nx,"ny":self.grid_ny,"dz":self.grid_dz,
                    "z_min":CACHE_Z_MIN,"z_max":CACHE_Z_MAX,"blend":CACHE_BLEND_WIDTH}
        identity = hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()
        self.cache_path = self.cache_dir/f"saip-{identity}.npz"
        padding = 3*self.grid_dz
        self._z_start = CACHE_Z_MIN-padding
        nz = int(np.ceil((CACHE_Z_MAX-CACHE_Z_MIN+2*padding)/self.grid_dz))+1
        expected_shape = (2, nz, self.grid_nx, self.grid_ny)
        if self.cache_path.is_file():
            try:
                with np.load(self.cache_path, allow_pickle=False) as stored:
                    grid = stored["potential"]
                    saved_settings = str(stored["settings"].item())
                    if grid.shape == expected_shape and np.isfinite(grid).all() and saved_settings == json.dumps(settings, sort_keys=True):
                        self._grid = grid
                        self.cache_loaded = True
                        return
            except (OSError, ValueError, KeyError, EOFError):
                # Rebuild a truncated/obsolete cache; never use partial data.
                pass
        start = time.perf_counter()
        grid = np.zeros(expected_shape)
        dx, dy = self.period/np.array([self.grid_nx, self.grid_ny])
        qx, qy = np.meshgrid(np.arange(-int(np.ceil(CUTOFF/dx)), int(np.ceil(CUTOFF/dx))+1),
                             np.arange(-int(np.ceil(CUTOFF/dy)), int(np.ceil(CUTOFF/dy))+1), indexing="ij")
        lateral_square = (qx*dx)**2+(qy*dy)**2
        flat_bins = (qx % self.grid_nx)*self.grid_ny+(qy % self.grid_ny)
        # Au basis and ABC layer shifts fall on grid nodes. Consequently one
        # finite-range lattice convolution yields all periodic query nodes,
        # rather than repeatedly rebuilding a giant gridpoint/Au pair array.
        for iz in range(nz):
            z = self._z_start+iz*self.grid_dz
            for layer, height in enumerate(self.layer_heights):
                dz = z-height
                keep = lateral_square+dz*dz < CUTOFF*CUTOFF
                if not np.any(keep):
                    continue
                radius = np.sqrt(lateral_square[keep]+dz*dz)
                if np.any(radius < 1e-8):
                    raise ValueError("A cache node coincides with an Au site; change the vertical cache spacing.")
                taper, _ = _taper(radius)
                bins = flat_bins[keep]
                for code, element in enumerate(("C", "H")):
                    p, gold = self.parameters[element, "Au"], self.parameters["Au", element]
                    damping = np.exp(-np.logaddexp(0, -p.d*(radius/p.damping_radius-1)))
                    kernel = .5*p.epsilon*np.exp(-p.decay*(radius-p.beta))
                    kernel += (.5*gold.epsilon+gold.c)*np.exp(-gold.decay*(radius-gold.beta))
                    kernel -= p.c6*radius**-6*damping
                    if element == "H":
                        kernel += p.c*np.exp(-p.decay*(radius-p.beta)-lateral_square[keep]/p.delta**2)
                    field = np.bincount(bins, weights=taper*kernel,
                                        minlength=self.grid_nx*self.grid_ny).reshape(self.grid_nx, self.grid_ny)
                    sx, sy = (layer % 3)*self.grid_nx//2, (layer % 3)*self.grid_ny//6
                    grid[code, iz] += np.roll(field, (sx, sy), axis=(0, 1))
                    grid[code, iz] += np.roll(field, (sx+self.grid_nx//2, sy+self.grid_ny//2), axis=(0, 1))
        self._grid = grid
        self.cache_build_seconds = time.perf_counter()-start
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        # The filename is fully determined by physics and interpolation data;
        # replacement is atomic even if two preparation workers share a cache.
        temporary = self.cache_path.with_name(self.cache_path.name+f".{os.getpid()}.{time.time_ns()}.tmp")
        try:
            with temporary.open("wb") as output:
                np.savez_compressed(output, potential=grid, settings=json.dumps(settings, sort_keys=True))
            os.replace(temporary, self.cache_path)
        finally:
            if temporary.exists():
                temporary.unlink()

    def _interpolate(self, positions, active):
        if not len(active):
            return np.empty(0), np.empty((0, 3))
        values = positions[active]
        spacing = np.array([*self.period/np.array([self.grid_nx, self.grid_ny]), self.grid_dz])
        scaled = values/spacing
        scaled[:, 2] = (values[:, 2]-self._z_start)/self.grid_dz
        base = np.floor(scaled).astype(np.int64)
        fraction = scaled-base
        wx, dx = _cubic_weights(fraction[:, 0])
        wy, dy = _cubic_weights(fraction[:, 1])
        wz, dz = _cubic_weights(fraction[:, 2])
        offset = np.arange(-1, 3)
        x = (base[:, 0, None]+offset) % self.grid_nx
        y = (base[:, 1, None]+offset) % self.grid_ny
        z = base[:, 2, None]+offset
        coefficients = self._grid[self._element_codes[active, None, None, None],
                                  z[:, :, None, None], x[:, None, :, None], y[:, None, None, :]]
        energy = np.einsum("nijk,ni,nj,nk->n", coefficients, wz, wx, wy)
        gradient = np.column_stack([
            np.einsum("nijk,ni,nj,nk->n", coefficients, wz, dx, wy),
            np.einsum("nijk,ni,nj,nk->n", coefficients, wz, wx, dy),
            np.einsum("nijk,ni,nj,nk->n", coefficients, dz, wx, wy)])
        return energy, -gradient/spacing

    def energy_and_forces_eV(self, positions, center_weights=None, center_dz=None):
        positions = self._positions(positions)
        if self._compiled_kernel is not None:
            return self._evaluate_compiled(positions, center_weights, center_dz)
        geometry = self.normal_geometry(positions)
        normals = geometry[0]
        count = len(positions)
        if self._grid is None:
            weight, weight_derivative = np.zeros(count), np.zeros(count)
        else:
            weight, weight_derivative = _cache_weight(positions[:, 2])
        cache_energy, cache_force = np.zeros(count), np.zeros_like(positions)
        cache_atoms = np.flatnonzero((weight > 0) & (self._element_codes >= 0))
        if len(cache_atoms):
            cache_energy[cache_atoms], cache_force[cache_atoms] = self._interpolate(positions, cache_atoms)
        # Each exact pair array is built once. Fixed-z H needs no pairs in the
        # interior cache domain; C and H with real local normals always retain
        # every direction-dependent pair up to the original 16 A cutoff.
        anisotropic_atoms = (self._element_codes == 0) | (geometry[2] >= 2)
        active = np.flatnonzero(((weight < 1) | anisotropic_atoms) & (self._element_codes >= 0))
        indices, delta, radius = self._pairs(positions, active)
        base_energy = np.zeros(count)
        base_force = np.zeros_like(positions)
        anis_energy = np.zeros(count)
        force = np.zeros_like(positions)
        normal_gradient = np.zeros_like(positions)
        for code, element in enumerate(("C", "H")):
            selected = self._element_codes[indices] == code
            if not np.any(selected):
                continue
            atom = indices[selected]
            displacement, distance = delta[selected], radius[selected]
            p, gold = self.parameters[element, "Au"], self.parameters["Au", element]
            exact = weight[atom] < 1
            if np.any(exact):
                exact_atoms = atom[exact]
                values, direct = _radial(displacement[exact], distance[exact], p, gold)
                if element == "H":
                    fixed = np.zeros_like(displacement[exact])
                    fixed[:, 2] = 1
                    fixed_energy, fixed_force, _ = _anisotropic(displacement[exact], distance[exact], fixed, p)
                    values += fixed_energy
                    direct += fixed_force
                base_energy += np.bincount(exact_atoms, weights=values, minlength=count)
                for axis in range(3):
                    base_force[:, axis] += np.bincount(exact_atoms, weights=direct[:, axis], minlength=count)
            directed = anisotropic_atoms[atom]
            if not np.any(directed):
                continue
            atom = atom[directed]
            displacement, distance = displacement[directed], distance[directed]
            anis, extra, gradient = _anisotropic(displacement, distance, normals[atom], p)
            if element == "H":
                fixed = np.zeros_like(displacement)
                fixed[:, 2] = 1
                fixed_energy, fixed_force, _ = _anisotropic(displacement, distance, fixed, p)
                anis -= fixed_energy
                extra -= fixed_force
            anis_energy += np.bincount(atom, weights=anis, minlength=count)
            for axis in range(3):
                force[:, axis] += np.bincount(atom, weights=extra[:, axis], minlength=count)
                normal_gradient[:, axis] += np.bincount(atom, weights=gradient[:, axis], minlength=count)
        force += (1-weight[:, None])*base_force+weight[:, None]*cache_force
        force[:, 2] -= weight_derivative*(cache_energy-base_energy)
        atomic_energy = (1-weight)*base_energy+weight*cache_energy+anis_energy
        if center_weights is not None:
            force *= center_weights[:, None]
            normal_gradient *= center_weights[:, None]
            force[:, 2] -= center_dz*atomic_energy
            atomic_energy *= center_weights
        _add_normal_chain_forces(force, normal_gradient, geometry)
        energy = np.sum(atomic_energy)
        return self._add_hetero((float(energy), force), positions, center_weights, center_dz)

    def energy_and_forces(self, positions):
        energy, force = self.energy_and_forces_eV(positions)
        return energy*EV_TO_KCAL_MOL, force*EV_TO_KCAL_MOL

    def metadata(self):
        result = {"model":"frozen six-layer periodic Au(111), SAIP CHAu",
            "algorithm_version":ALGORITHM_VERSION,"evaluation":"periodic exact lattice sum" if self._grid is None else "cached radial/H field + exact full-cutoff normal term",
            "parameter_file":str(self.parameter_path.resolve()),"parameter_sha256":self.parameter_sha256,
            "citation":"https://doi.org/10.1021/acs.jctc.1c00622",
            "supported_elements":["C","H"],"units":{"energy":"kcal/mol","force":"kcal/mol/angstrom","length":"angstrom"},
            "eV_API":"energy_and_forces_eV","cutoff_angstrom":CUTOFF,"taper":True,
            "surface_period_xy_angstrom":self.period.tolist(),"layer_heights_angstrom":self.layer_heights.tolist(),
            "reference_source":self.reference_source,"normal_policy":self.normal_policy,
            "normal_neighbours":("immutable covalent topology and preparation-time plane/surface-axis classes; regularized conservative chain forces"
                if self.normal_policy == BONDED_NORMAL_POLICY else "same-molecule unwrapped distances, at most 3; full normal chain forces"),
            "normal_model_is_game_adaptation":self.normal_policy == BONDED_NORMAL_POLICY,
            "normal_topology_sha256":self.normal_topology_sha256,
            "normal_modes":list(self.normal_modes) if self.normal_modes is not None else None,
            "normal_area_regularization_angstrom2":NORMAL_AREA_EPSILON if self.normal_policy == BONDED_NORMAL_POLICY else 0.,
            "candidate_Au_sites":len(self._tile),"cache_requested":self.use_cache,"cache_loaded":self.cache_loaded,
            "cache_build_seconds":self.cache_build_seconds,"preparation_seconds":self.preparation_seconds,
            "cache_file":str(self.cache_path) if self.cache_path else None,
            "cache_grid_shape":list(self._grid.shape) if self._grid is not None else None,
            "cache_spacing_angstrom":[self.period[0]/self.grid_nx,self.period[1]/self.grid_ny,self.grid_dz] if self._grid is not None else None,
            "cache_height_support_angstrom":[CACHE_Z_MIN,CACHE_Z_MAX] if self._grid is not None else None,
            "cache_blend_width_angstrom":CACHE_BLEND_WIDTH if self._grid is not None else None,
            "cache_interpolation":"Catmull-Rom scalar potential with analytic gradient; smooth conservative exact fallback" if self._grid is not None else None,
            "accelerator_requested":self.accelerator_requested,"accelerator":self.accelerator,
            "numba_version":self.numba_version,"numba_import_seconds":self.numba_import_seconds,
            "numba_compile_warmup_seconds":self.numba_compile_warmup_seconds,
            "compiled_execution":{"threads":1,"fastmath":False,"disk_cache":True} if self.accelerator == "numba" else None,
            "limitations":["Frozen laterally perfect Au; no phonons, reconstruction, steps or electrons.",
                ("Bonded sp/sp3 normal handling is a game adaptation, not validated SAIP parameter transferability for arbitrary C/H molecules."
                 if self.normal_policy == BONDED_NORMAL_POLICY else "C/H geometry and parameter transferability follow SAIP limitations; other elements and >3 normal neighbours rejected."),
                "Reference Au layer heights retained; tiny in-plane/nonuniform relaxation is replaced by a perfect periodic lattice."]}
        if self.hetero_surface is not None:
            result.update(model="frozen Au(111), CHAu SAIP + provisional X-Au game LJ",
                          supported_elements=list(SUPPORTED_RACE_ELEMENTS),
                          hetero_au=self.hetero_surface.metadata())
            result["limitations"].extend(self.hetero_surface.parameters["limitations"])
        return result


def create_periodic_surface(atoms, reference_positions=None, **kwargs):
    return PeriodicAuPotential(atoms, reference_positions, **kwargs)
