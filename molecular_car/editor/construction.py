"""Initial editor geometry, deliberately independent of force fields.

Existing conformers and the native C60 template are boundaries, not molecules
to re-embed. Only a genuinely new connected fragment is embedded, with one
local cap atom if it has a single attachment. Array indices carry identity.
"""
from __future__ import annotations

import math

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D


class ConstructionError(ValueError):
    pass


def embed_initial(molecule, random_seed):
    """ETKDG initial coordinates, with no UFF requirement or minimization."""
    params = AllChem.ETKDGv3()
    params.randomSeed, params.numThreads = random_seed, 1
    params.maxIterations, params.timeout = 200, 10
    with rdBase.BlockLogs():
        status = AllChem.EmbedMolecule(molecule, params)
    if status != 0:
        params.useRandomCoords = True
        with rdBase.BlockLogs():
            status = AllChem.EmbedMolecule(molecule, params)
    if status != 0:
        raise ConstructionError("无法生成符合连接及立体标记的初始三维构象；请检查新增结构。")
    return molecule


def _unit(vector, fallback=(1.0, 0.0, 0.0)):
    vector = np.asarray(vector, dtype=float)
    length = float(np.linalg.norm(vector))
    return vector / length if length > 1e-10 else np.asarray(fallback, dtype=float)


def _perpendicular(vector):
    axis = np.eye(3)[np.argmin(np.abs(vector))]
    return _unit(np.cross(vector, axis))


def _rotation(source, target):
    """A proper rotation; fragment alignment must never reflect chirality."""
    source, target = _unit(source), _unit(target)
    cosine = float(np.clip(np.dot(source, target), -1, 1))
    cross = np.cross(source, target)
    sine = float(np.linalg.norm(cross))
    if sine < 1e-10:
        if cosine > 0:
            return np.eye(3)
        axis = _perpendicular(source)
        return 2 * np.outer(axis, axis) - np.eye(3)
    x, y, z = cross
    skew = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    return np.eye(3) + skew + skew @ skew * ((1 - cosine) / sine ** 2)


def _axis_rotation(axis, angle):
    axis = _unit(axis)
    x, y, z = axis
    skew = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    return np.eye(3) + math.sin(angle) * skew + (1 - math.cos(angle)) * (skew @ skew)


def _bond_length(molecule, first, second):
    table = Chem.GetPeriodicTable()
    radius = sum(table.GetRcovalent(molecule.GetAtomWithIdx(i).GetAtomicNum()) for i in (first, second))
    order = molecule.GetBondBetweenAtoms(first, second).GetBondTypeAsDouble()
    factor = {1: 1.0, 1.5: 0.92, 2: 0.9, 3: 0.8}.get(order, 1.0)
    return max(0.7, radius * factor)


def _attachment_direction(molecule, anchor, excluded, cages):
    points = molecule.GetConformer().GetPositions()
    for members in cages:
        if anchor in members:
            return _unit(points[anchor] - points[list(members)].mean(axis=0))
    vectors = [_unit(points[neighbor.GetIdx()] - points[anchor])
               for neighbor in molecule.GetAtomWithIdx(anchor).GetNeighbors()
               if neighbor.GetIdx() not in excluded
               and np.linalg.norm(points[neighbor.GetIdx()] - points[anchor]) > 0.05]
    if not vectors:
        return np.array([1.0, 0.0, 0.0])
    atom = molecule.GetAtomWithIdx(anchor)
    planar = atom.GetIsAromatic() or atom.GetHybridization() == Chem.HybridizationType.SP2
    if atom.GetHybridization() == Chem.HybridizationType.SP:
        return -vectors[0]
    if len(vectors) == 1:
        cosine = -0.5 if planar else -1 / 3
        return cosine * vectors[0] + math.sqrt(1 - cosine ** 2) * _perpendicular(vectors[0])
    if len(vectors) == 2 and not planar:
        # One of the two unoccupied tetrahedral directions, never a flat
        # continuation of provisional drawing coordinates.
        bisector = _unit(vectors[0] + vectors[1], _perpendicular(vectors[0]))
        normal = _unit(np.cross(vectors[0], vectors[1]), _perpendicular(bisector))
        cosine = min(0.99, (1 / 3) / max(0.34, float(np.dot(bisector, vectors[0]))))
        return _unit(-cosine * bisector + math.sqrt(1 - cosine ** 2) * normal)
    return _unit(-np.sum(vectors, axis=0), _perpendicular(vectors[0]))


def _components(molecule, indices):
    remaining = set(indices)
    while remaining:
        component = {min(remaining)}
        frontier = list(component)
        for index in frontier:
            for neighbor in molecule.GetAtomWithIdx(index).GetNeighbors():
                second = neighbor.GetIdx()
                if second in remaining and second not in component:
                    component.add(second)
                    frontier.append(second)
        remaining.difference_update(component)
        yield component


def _boundary(molecule, component):
    return [(index, neighbor.GetIdx()) for index in sorted(component)
            for neighbor in molecule.GetAtomWithIdx(index).GetNeighbors()
            if neighbor.GetIdx() not in component]


def _fragment(molecule, component, anchor):
    """Extract by exact indices, retaining every new center's neighbors."""
    included = sorted(component | ({anchor} if anchor is not None else set()))
    local_by_global = {index: local for local, index in enumerate(included)}
    builder = Chem.RWMol(molecule)
    for index in reversed(range(molecule.GetNumAtoms())):
        if index not in local_by_global:
            builder.RemoveAtom(index)
    if anchor is not None:
        # Remote atoms are not copied. A neutral cap supplies an attachment
        # vector without imposing the parent's ring or force-field parameters.
        cap = Chem.Atom(molecule.GetAtomWithIdx(anchor).GetAtomicNum())
        builder.ReplaceAtom(local_by_global[anchor], cap)
    fragment = builder.GetMol()
    fragment.RemoveAllConformers()
    with rdBase.BlockLogs():
        try:
            Chem.SanitizeMol(fragment)
        except Exception as exc:
            raise ConstructionError("新增片段不能独立生成初始构象；复杂共环或多点连接请手工构造。") from exc
    for index in component:
        original = molecule.GetAtomWithIdx(index)
        if original.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            if any(neighbor.GetIdx() not in local_by_global for neighbor in original.GetNeighbors()):
                raise ConstructionError("新增立体中心的连接信息不完整。请补齐连接或手动调整构象。")
    for bond in molecule.GetBonds():
        if (bond.GetBeginAtomIdx() in component or bond.GetEndAtomIdx() in component) and bond.GetStereo() not in (
                Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY):
            if any(index not in local_by_global for index in bond.GetStereoAtoms()):
                raise ConstructionError("新增立体键需要多个原有连接点；请手动调整构象或重建三维构象。")
    with rdBase.BlockLogs():
        fragment = Chem.AddHs(fragment)
    return fragment, local_by_global


def _choose_roll(molecule, component, anchor, endpoint, points, direction):
    """Choose among rigid rotations using clearance, not an optimizer."""
    original = molecule.GetConformer().GetPositions()
    old = [i for i in range(molecule.GetNumAtoms()) if i not in component
           and molecule.GetAtomWithIdx(i).GetAtomicNum() != 1]
    heavy = [i for i in sorted(component) if molecule.GetAtomWithIdx(i).GetAtomicNum() != 1]
    if not old or not heavy:
        return points
    pivot = points[endpoint]
    best, best_score = points, -math.inf
    for step in range(12):
        rotation = _axis_rotation(direction, step * math.pi / 6)
        candidate = {index: (point - pivot) @ rotation.T + pivot for index, point in points.items()}
        scores = [np.linalg.norm(candidate[i] - original[j])
                  for i in heavy for j in old if not (i == endpoint and j == anchor)]
        score = min(scores) if scores else math.inf
        if score > best_score:
            best, best_score = candidate, score
    return best


def construct_new_regions(molecule, new_indices, seed_groups, cages, fixed_indices, random_seed):
    """Place rigid cages and locally build genuinely new single-point fragments.

    The input conformer is updated only at new, unpinned indices. Existing
    atoms (including manual and previously generated H) are never moved.
    Complex bridges fail before the caller applies this temporary result.
    """
    conformer = molecule.GetConformer()
    new_indices = set(new_indices) - set(fixed_indices)
    seeded = {index for group in seed_groups for index in group}
    for group in seed_groups:
        group = set(group)
        external = [(inside, outside) for inside, outside in _boundary(molecule, group)
                    if outside not in new_indices]
        if len(external) > 1:
            raise ConstructionError("C60 新增球壳具有多个原有固定连接点；无法自动放置。请手动调整构象或运行 UFF 优化。")
        if external and not group.issubset(new_indices):
            raise ConstructionError("C60 新增球壳的固定原子与外部连接同时限制放置；请手工调整连接。")
        if external:
            endpoint, anchor = external[0]
            original = conformer.GetPositions()
            outward = _attachment_direction(molecule, anchor, new_indices, cages)
            radial = original[endpoint] - original[list(group)].mean(axis=0)
            rotation = _rotation(radial, -outward)
            target = original[anchor] + _bond_length(molecule, endpoint, anchor) * outward
            for index in group:
                point = (original[index] - original[endpoint]) @ rotation.T + target
                conformer.SetAtomPosition(index, Point3D(*point))
    built = []
    for component in _components(molecule, new_indices - seeded):
        boundary = _boundary(molecule, component)
        if len(boundary) > 1:
            raise ConstructionError("新增片段具有多个原有固定连接点；无法自动放置。请手动调整构象或重建三维构象。")
        endpoint, anchor = boundary[0] if boundary else (None, None)
        original = conformer.GetPositions()
        if len(component) == 1:
            index = next(iter(component))
            if anchor is not None:
                direction = _attachment_direction(molecule, anchor, component, cages)
                point = original[anchor] + _bond_length(molecule, index, anchor) * direction
                conformer.SetAtomPosition(index, Point3D(*point))
            built.append({"atomCount": 1, "method": "local-valence-geometry"})
            continue
        fragment, local = _fragment(molecule, component, anchor)
        embed_initial(fragment, random_seed)
        embedded = fragment.GetConformer().GetPositions()
        if anchor is not None:
            direction = _attachment_direction(molecule, anchor, component, cages)
            rotation = _rotation(embedded[local[endpoint]] - embedded[local[anchor]], direction)
            target = original[anchor] + _bond_length(molecule, endpoint, anchor) * direction
            points = {index: (embedded[local[index]] - embedded[local[endpoint]]) @ rotation.T + target
                      for index in component}
            points = _choose_roll(molecule, component, anchor, endpoint, points, direction)
        else:
            translation = original[list(component)].mean(axis=0) - embedded[[local[i] for i in component]].mean(axis=0)
            points = {index: embedded[local[index]] + translation for index in component}
        for index, point in points.items():
            conformer.SetAtomPosition(index, Point3D(*point))
        built.append({"atomCount": len(component), "method": "local-ETKDG"})
    if not np.isfinite(conformer.GetPositions()).all():
        raise ConstructionError("局部初始构象含无效坐标；")
    return built
