"""Single-thread compiled evaluation of the selected periodic SAIP policy.

The scalar loops remove temporary pair arrays and Python dispatch overhead.
It follows either the distance-reference or explicitly bonded-game normal
policy, retaining its directional term and complete analytic force derivative.
The radial interpolation has the same conservative fallback as periodic_au.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(__file__).resolve().parents[1]/".local/numba-cache"))
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import numba
import numpy as np

NUMBA_VERSION = numba.__version__


@numba.njit(cache=True, fastmath=False, inline="always")
def _taper(radius):
    x = radius/16.0
    square = x*x
    cube = square*x
    offset = x-1.0
    return (1+square*square*(-35+x*(84+x*(-70+20*x))),
            140*cube*offset*offset*offset/16.0)


@numba.njit(cache=True, fastmath=False, inline="always")
def _radial(radius, p, gold, taper, taper_derivative):
    # Parameter columns: beta, decay, delta^2, epsilon, C, d, sr*reff, C6.
    organic_rep = .5*p[3]*math.exp(-p[1]*(radius-p[0]))
    gold_rep = (.5*gold[3]+gold[4])*math.exp(-gold[1]*(radius-gold[0]))
    argument = p[5]*(radius/p[6]-1)
    if argument >= 0:
        damping = 1/(1+math.exp(-argument))
    else:
        exponential = math.exp(argument)
        damping = exponential/(1+exponential)
    dispersion = -p[7]*radius**-6*damping
    raw = organic_rep+gold_rep+dispersion
    derivative = -p[1]*organic_rep-gold[1]*gold_rep
    derivative += p[7]*radius**-7*damping*(6-radius*p[5]/p[6]*(1-damping))
    return taper*raw, -(taper_derivative*raw+taper*derivative)/radius


@numba.njit(cache=True, fastmath=False, inline="always")
def _anisotropic(dx, dy, dz, radius, nx, ny, nz, p, taper, taper_derivative):
    projection = dx*nx+dy*ny+dz*nz
    px, py, pz = dx-projection*nx, dy-projection*ny, dz-projection*nz
    rho_square = dx*dx+dy*dy+dz*dz-projection*projection
    exponential = p[4]*math.exp(-p[1]*(radius-p[0])-rho_square/p[2])
    energy = taper*exponential
    radial_force = -(taper_derivative-p[1]*taper)*exponential/radius
    lateral_force = 2*energy/p[2]
    normal_force = lateral_force*projection
    return (energy, radial_force*dx+lateral_force*px,
            radial_force*dy+lateral_force*py, radial_force*dz+lateral_force*pz,
            normal_force*dx, normal_force*dy, normal_force*dz)


@numba.njit(cache=True, fastmath=False, inline="always")
def _weights(t):
    square, cube = t*t, t*t*t
    return ((-.5*t+square-.5*cube, 1-2.5*square+1.5*cube,
             .5*t+2*square-1.5*cube, -.5*square+.5*cube),
            (-.5+2*t-1.5*square, -5*t+4.5*square,
             .5+4*t-4.5*square, -t+1.5*square))


@numba.njit(cache=True, fastmath=False, inline="always")
def _interpolate(grid, code, x, y, z, period_x, period_y, grid_dz, z_start):
    spacing_x, spacing_y = period_x/grid.shape[2], period_y/grid.shape[3]
    sx, sy, sz = (x % period_x)/spacing_x, (y % period_y)/spacing_y, (z-z_start)/grid_dz
    bx, by, bz = math.floor(sx), math.floor(sy), math.floor(sz)
    wx, dwx = _weights(sx-bx)
    wy, dwy = _weights(sy-by)
    wz, dwz = _weights(sz-bz)
    energy, gx, gy, gz = 0.0, 0.0, 0.0, 0.0
    for k in range(4):
        zi = bz+k-1
        for i in range(4):
            xi = (bx+i-1) % grid.shape[2]
            for j in range(4):
                yi = (by+j-1) % grid.shape[3]
                value = grid[code, zi, xi, yi]
                energy += value*wz[k]*wx[i]*wy[j]
                gx += value*wz[k]*dwx[i]*wy[j]
                gy += value*wz[k]*wx[i]*dwy[j]
                gz += value*dwz[k]*wx[i]*wy[j]
    return energy, -gx/spacing_x, -gy/spacing_y, -gz/grid_dz


@numba.njit(cache=True, fastmath=False, inline="always")
def _cache_weight(z):
    low = min(1.0, max(0.0, (z-2.0)/.3))
    high = min(1.0, max(0.0, (12.0-z)/.3))
    low_weight = low**3*(10+low*(-15+6*low))
    high_weight = high**3*(10+high*(-15+6*high))
    low_derivative = 30*low*low*(low-1)**2/.3
    high_derivative = 30*high*high*(high-1)**2/.3
    return low_weight*high_weight, low_derivative*high_weight-low_weight*high_derivative


@numba.njit(cache=True, fastmath=False)
def _normal_geometry(positions, normal_cutoffs, topology_neighbours, topology_counts, normal_epsilon):
    count = len(positions)
    normals = np.zeros((count, 3))
    neighbours = np.full((count, 3), -1, dtype=np.int64)
    neighbour_count = np.zeros(count, dtype=np.int64)
    first = np.zeros((count, 3))
    second = np.zeros((count, 3))
    normal_length = np.ones(count)
    for i in range(count):
        normals[i, 2] = 1.0
        if normal_epsilon > 0:
            neighbour_count[i] = topology_counts[i]
            for k in range(3):
                neighbours[i, k] = topology_neighbours[i, k]
        else:
            for j in range(count):
                dx, dy, dz = (positions[j, 0]-positions[i, 0],
                              positions[j, 1]-positions[i, 1], positions[j, 2]-positions[i, 2])
                square = dx*dx+dy*dy+dz*dz
                cutoff = normal_cutoffs[i, j]
                if square > 0 and square < cutoff*cutoff:
                    n = neighbour_count[i]
                    if n == 3:
                        raise ValueError("A SAIP atom has more than three normal neighbours; this geometry is unsupported.")
                    neighbours[i, n] = j
                    neighbour_count[i] = n+1
        if neighbour_count[i] >= 2:
            pivot = neighbours[i, 2] if neighbour_count[i] == 3 else i
            for axis in range(3):
                first[i, axis] = positions[neighbours[i, 0], axis]-positions[pivot, axis]
                second[i, axis] = positions[neighbours[i, 1], axis]-positions[pivot, axis]
            ax, ay, az = first[i, 0], first[i, 1], first[i, 2]
            bx, by, bz = second[i, 0], second[i, 1], second[i, 2]
            rx, ry, rz = ay*bz-az*by, az*bx-ax*bz, ax*by-ay*bx
            length = math.sqrt(rx*rx+ry*ry+rz*rz+normal_epsilon*normal_epsilon)
            if normal_epsilon == 0 and length < 1e-12:
                raise ValueError("A SAIP local plane has zero normal magnitude.")
            normal_length[i] = length
            normals[i, 0], normals[i, 1], normals[i, 2] = rx/length, ry/length, rz/length
    return normals, neighbours, neighbour_count, first, second, normal_length


@numba.njit(cache=True, fastmath=False)
def energy_and_forces(positions, codes, normal_cutoffs, tile, parameters, grid,
                      period_x, period_y, grid_dz, z_start, has_cache,
                      topology_neighbours, topology_counts, normal_epsilon,
                      center_weights=None, center_dz=None):
    """Return eV and eV/A for complete unwrapped coordinates.

    Parameter data and cached grid are independent of molecule coordinates.
    This function has no Python callbacks, parallel loops, or fastmath flags.
    """
    count = len(positions)
    normals, neighbours, neighbour_count, first, second, normal_length = _normal_geometry(
        positions, normal_cutoffs, topology_neighbours, topology_counts, normal_epsilon)
    normal_gradient = np.zeros((count, 3))
    force = np.zeros((count, 3))
    total_energy = 0.0
    for i in range(count):
        code = codes[i]
        if code < 0:
            continue  # independent isotropic X-Au potential; keep full normal graph
        p, gold = parameters[code, 0], parameters[code, 1]
        weight, derivative = _cache_weight(positions[i, 2]) if has_cache else (0.0, 0.0)
        cached_e, cached_fx, cached_fy, cached_fz = 0.0, 0.0, 0.0, 0.0
        if weight > 0:
            cached_e, cached_fx, cached_fy, cached_fz = _interpolate(
                grid, code, positions[i, 0], positions[i, 1], positions[i, 2],
                period_x, period_y, grid_dz, z_start)
        base_e, base_fx, base_fy, base_fz = 0.0, 0.0, 0.0, 0.0
        anis_e, anis_fx, anis_fy, anis_fz = 0.0, 0.0, 0.0, 0.0
        direction_dependent = code == 0 or neighbour_count[i] >= 2
        if weight < 1 or direction_dependent:
            x, y, z = positions[i, 0] % period_x, positions[i, 1] % period_y, positions[i, 2]
            for site in range(len(tile)):
                dx, dy, dz = x-tile[site, 0], y-tile[site, 1], z-tile[site, 2]
                square = dx*dx+dy*dy+dz*dz
                if square >= 256.0:
                    continue
                if square < 1e-16:
                    raise ValueError("A molecular atom coincides with a gold lattice site.")
                radius = math.sqrt(square)
                taper, taper_derivative = _taper(radius)
                fixed_e, fixed_fx, fixed_fy, fixed_fz = 0.0, 0.0, 0.0, 0.0
                if code == 1:
                    fixed = _anisotropic(dx, dy, dz, radius, 0.0, 0.0, 1.0, p, taper, taper_derivative)
                    fixed_e, fixed_fx, fixed_fy, fixed_fz = fixed[0], fixed[1], fixed[2], fixed[3]
                if weight < 1:
                    radial_e, radial_force = _radial(radius, p, gold, taper, taper_derivative)
                    base_e += radial_e+fixed_e
                    base_fx += radial_force*dx+fixed_fx
                    base_fy += radial_force*dy+fixed_fy
                    base_fz += radial_force*dz+fixed_fz
                if direction_dependent:
                    actual = _anisotropic(dx, dy, dz, radius, normals[i, 0], normals[i, 1], normals[i, 2],
                                          p, taper, taper_derivative)
                    anis_e += actual[0]-fixed_e
                    anis_fx += actual[1]-fixed_fx
                    anis_fy += actual[2]-fixed_fy
                    anis_fz += actual[3]-fixed_fz
                    normal_gradient[i, 0] += actual[4]
                    normal_gradient[i, 1] += actual[5]
                    normal_gradient[i, 2] += actual[6]
        atom_energy = (1-weight)*base_e+weight*cached_e+anis_e
        total_energy += atom_energy if center_weights is None else center_weights[i]*atom_energy
        force[i, 0] += (1-weight)*base_fx+weight*cached_fx+anis_fx
        force[i, 1] += (1-weight)*base_fy+weight*cached_fy+anis_fy
        force[i, 2] += (1-weight)*base_fz+weight*cached_fz+anis_fz-derivative*(cached_e-base_e)
        if center_weights is not None:
            for axis in range(3):
                force[i, axis] *= center_weights[i]
                normal_gradient[i, axis] *= center_weights[i]
            force[i, 2] -= center_dz[i]*atom_energy
    # Reverse-mode derivative of the normalized local plane. The molecule
    # is never folded into an Au cell while finding or applying this chain.
    for i in range(count):
        if neighbour_count[i] < 2:
            continue
        nx, ny, nz = normals[i, 0], normals[i, 1], normals[i, 2]
        dot = normal_gradient[i, 0]*nx+normal_gradient[i, 1]*ny+normal_gradient[i, 2]*nz
        gx = (normal_gradient[i, 0]-dot*nx)/normal_length[i]
        gy = (normal_gradient[i, 1]-dot*ny)/normal_length[i]
        gz = (normal_gradient[i, 2]-dot*nz)/normal_length[i]
        ax, ay, az = first[i, 0], first[i, 1], first[i, 2]
        bx, by, bz = second[i, 0], second[i, 1], second[i, 2]
        gax, gay, gaz = by*gz-bz*gy, bz*gx-bx*gz, bx*gy-by*gx
        gbx, gby, gbz = gy*az-gz*ay, gz*ax-gx*az, gx*ay-gy*ax
        n0, n1 = neighbours[i, 0], neighbours[i, 1]
        pivot = neighbours[i, 2] if neighbour_count[i] == 3 else i
        force[n0, 0] -= gax
        force[n0, 1] -= gay
        force[n0, 2] -= gaz
        force[n1, 0] -= gbx
        force[n1, 1] -= gby
        force[n1, 2] -= gbz
        force[pivot, 0] += gax+gbx
        force[pivot, 1] += gay+gby
        force[pivot, 2] += gaz+gbz
    return total_energy, force
