"""Independent Au(111) reference experiment, never used by the game engine.

Coordinates use angstroms. The ABC-stacked rectangular surface cell is
constructed analytically from FCC geometry. Its default lattice constant is
the value in LAMMPS's bundled Au_u3.eam (Foiles, Baskes and Daw, PRB 33,
7983, 1986, https://doi.org/10.1103/PhysRevB.33.7983).

This module supplies geometry and output helpers, not a carbon/gold force
field. EAM applies only to gold. Cross interactions must be specified
independently and their chemical scope documented by the calling experiment.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import hashlib

import numpy as np


@dataclass
class AuSlab:
    positions: np.ndarray
    layers: np.ndarray
    box: np.ndarray
    lattice_constant: float
    nx: int
    ny: int
    layer_count: int
    frozen_layers: int

    @property
    def frozen_count(self) -> int:
        return int(np.count_nonzero(self.layers < self.frozen_layers))

    @property
    def nearest_spacing(self) -> float:
        return self.lattice_constant / np.sqrt(2.0)

    @property
    def lengths_xy(self) -> np.ndarray:
        return self.box[:2, 1] - self.box[:2, 0]

    def with_positions(self, positions: np.ndarray, *, top_zero: bool = False) -> AuSlab:
        coordinates = np.asarray(positions, dtype=float).copy()
        if coordinates.shape != self.positions.shape or not np.isfinite(coordinates).all():
            raise ValueError("Gold coordinates must be finite and preserve the atom count.")
        if top_zero:
            coordinates[:, 2] -= np.mean(coordinates[self.layers == self.layer_count - 1, 2])
        return replace(self, positions=coordinates)


def build_au111(nx: int = 20, ny: int = 12, layers: int = 6,
                frozen_layers: int = 2, lattice_constant: float = 4.08) -> AuSlab:
    """Make a periodic orthogonal FCC(111) slab, top layer at z=0.

    Each rectangular cell contains two atoms per layer. Surface vectors are
    [1 -1 0] and [1 1 -2]; z is [1 1 1]. Layer labels and atom IDs are ordered
    bottom to top, so the fixed bottom layers form a contiguous ID interval.
    """
    for value, name in ((nx, "nx"), (ny, "ny"), (layers, "layers"),
                        (frozen_layers, "frozen_layers")):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer.")
    if frozen_layers >= layers:
        raise ValueError("At least one gold layer must remain mobile.")
    if not np.isfinite(lattice_constant) or lattice_constant <= 0:
        raise ValueError("The lattice constant must be finite and positive.")
    d = lattice_constant / np.sqrt(2.0)
    h = np.sqrt(3.0) * d
    gap = lattice_constant / np.sqrt(3.0)
    lengths = np.array([nx * d, ny * h])
    coordinates = []
    labels = []
    for layer in range(layers):
        shift = (layer % 3) * np.array([d / 2.0, h / 6.0])
        for row in range(ny):
            for column in range(nx):
                for basis in ((0.0, 0.0), (d / 2.0, h / 2.0)):
                    xy = np.mod(np.array([column * d, row * h]) + basis + shift, lengths)
                    coordinates.append([*xy, (layer - layers + 1) * gap])
                    labels.append(layer)
    box = np.array([[0.0, lengths[0]], [0.0, lengths[1]],
                    [-(layers - 1) * gap - 20.0, 30.0]])
    return AuSlab(np.array(coordinates), np.array(labels, dtype=int), box,
                  float(lattice_constant), nx, ny, layers, frozen_layers)


def minimum_image_xy(displacements: np.ndarray, lengths_xy: np.ndarray) -> np.ndarray:
    result = np.asarray(displacements, dtype=float).copy()
    result[..., :2] -= np.rint(result[..., :2] / lengths_xy) * lengths_xy
    return result


def layer_statistics(slab: AuSlab) -> dict:
    rows = []
    for layer in range(slab.layer_count):
        z = slab.positions[slab.layers == layer, 2]
        rows.append({"layer_bottom_to_top": layer, "count": len(z),
                     "mean_z_angstrom": float(np.mean(z)),
                     "std_z_angstrom": float(np.std(z)),
                     "min_z_angstrom": float(np.min(z)),
                     "max_z_angstrom": float(np.max(z)),
                     "frozen": layer < slab.frozen_layers})
    means = [row["mean_z_angstrom"] for row in rows]
    return {"layers": rows, "interlayer_gaps_angstrom": np.diff(means).tolist()}


def potential_metadata(path: Path) -> dict:
    payload = path.read_bytes()
    return {"path": str(path.resolve()), "filename": path.name,
            "sha256": hashlib.sha256(payload).hexdigest(),
            "header": payload.decode("utf-8", errors="replace").splitlines()[:3]}


def write_atomic_data(path: Path, slab: AuSlab) -> None:
    lines = ["Analytical ABC Au(111); units metal; atom_style atomic", "",
             f"{len(slab.positions)} atoms", "1 atom types", ""]
    for bounds, names in zip(slab.box, ("xlo xhi", "ylo yhi", "zlo zhi")):
        lines.append(f"{bounds[0]:.14g} {bounds[1]:.14g} {names}")
    lines.extend(["", "Masses", "", "1 196.97", "", "Atoms # atomic", ""])
    lines.extend(f"{index} 1 {x:.14g} {y:.14g} {z:.14g}"
                 for index, (x, y, z) in enumerate(slab.positions, 1))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_xyz(path: Path, slab: AuSlab) -> None:
    lines = [str(len(slab.positions)), "Au(111), angstrom, periodic x/y; top mean z=0"]
    lines.extend(f"Au {x:.10f} {y:.10f} {z:.10f}" for x, y, z in slab.positions)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_height_svg(path: Path, slab: AuSlab) -> None:
    """Standalone top-layer map; point colors encode actual height differences."""
    top = slab.positions[slab.layers == slab.layer_count - 1]
    width, height = 700.0, 700.0 * slab.lengths_xy[1] / slab.lengths_xy[0]
    deviations = top[:, 2] - np.mean(top[:, 2])
    scale = max(float(np.max(np.abs(deviations))), 0.001)
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height + 70}">',
             '<rect width="100%" height="100%" fill="#101522"/>',
             '<text x="20" y="24" fill="white" font-size="18">Relaxed Au(111) top layer</text>',
             f'<text x="20" y="48" fill="white" font-size="13">Height deviation: '
             f'{float(np.min(deviations)):.6g} to {float(np.max(deviations)):.6g} Å</text>']
    for (x, y, _), delta in zip(top, deviations):
        signed = delta / scale
        rgb = (int(190 + 55 * signed), int(160 - 40 * abs(signed)), int(140 - 70 * signed))
        lines.append(f'<circle cx="{x / slab.lengths_xy[0] * width:.3f}" '
                     f'cy="{70 + (1 - y / slab.lengths_xy[1]) * height:.3f}" r="5" '
                     f'fill="rgb{rgb}"><title>z deviation {delta:.9g} Å</title></circle>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")
