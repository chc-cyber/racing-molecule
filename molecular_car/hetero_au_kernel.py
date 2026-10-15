"""Compiled conservative scalar UFF-vdW/Au lattice extension, in eV and A."""
import math
import numpy as np
from .periodic_au_kernel import numba, _taper, _cache_weight, _interpolate


@numba.njit(cache=True, fastmath=False)
def energy_and_forces(positions, indices, coefficients, tile, grid, period_x,
                      period_y, grid_dz, z_start, has_cache, center_weights=None, center_dz=None):
    forces = np.zeros_like(positions)
    energy = 0.
    for k in range(len(indices)):
        i = indices[k]
        x, y, z = positions[i, 0] % period_x, positions[i, 1] % period_y, positions[i, 2]
        a, b = coefficients[k, 0], coefficients[k, 1]
        w, dw = _cache_weight(z) if has_cache else (0., 0.)
        ce, cf = 0., np.zeros(3)
        if w > 0:
            for code in range(2):
                scale = a if code == 0 else b
                e, fx, fy, fz = _interpolate(grid, code, x, y, z, period_x, period_y, grid_dz, z_start)
                ce += scale*e
                cf[0] += scale*fx; cf[1] += scale*fy; cf[2] += scale*fz
        ee, ef = 0., np.zeros(3)
        if w < 1:
            for j in range(len(tile)):
                dx, dy, dz = x-tile[j, 0], y-tile[j, 1], z-tile[j, 2]
                r2 = dx*dx+dy*dy+dz*dz
                if r2 >= 256.:
                    continue
                if r2 < 1e-16:
                    raise ValueError("A molecular atom coincides with a gold lattice site.")
                r = math.sqrt(r2)
                t, dt = _taper(r)
                r6 = r2**-3
                raw = a*r6*r6+b*r6
                derivative = (-12*a*r6*r6-6*b*r6)/r
                factor = -(dt*raw+t*derivative)/r
                ee += t*raw
                ef[0] += factor*dx; ef[1] += factor*dy; ef[2] += factor*dz
        atom_energy = w*ce+(1-w)*ee
        energy += atom_energy if center_weights is None else center_weights[i]*atom_energy
        for axis in range(3):
            forces[i, axis] = w*cf[axis]+(1-w)*ef[axis]
        forces[i, 2] -= dw*(ce-ee)
        if center_weights is not None:
            for axis in range(3):forces[i, axis] *= center_weights[i]
            forces[i, 2] -= center_dz[i]*atom_energy
    return energy, forces
