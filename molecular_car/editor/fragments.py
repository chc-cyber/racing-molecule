"""Procedurally generated chemical fragments, not downloaded visual meshes.

C60 is the graph of a truncated icosahedron: each oriented icosahedron edge
becomes an atom; reversed edge pairs are its thirty 6:6 double bonds. Adjacent
edges at each vertex form twelve single-bond pentagons. The resulting graph has
60 vertices, 90 edges, degree three and bond-order sum four at every carbon.
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem

from .service import FORMAT, _layout


def _base(name, slug, points, edges, max_order=1, rigid=False):
    atoms = [{"id": f"{slug}-a{i + 1}", "element": "C", "charge": 0,
              "chiralTag": "CHI_UNSPECIFIED", "pos2": list(point[:2]), "pos3": list(point)}
             for i, point in enumerate(points)]
    bonds = [{"id": f"{slug}-b{i + 1}", "a": atoms[a]["id"], "b": atoms[b]["id"],
              "order": order, "stereo": "STEREONONE", "direction": "NONE", "stereoAtoms": []}
             for i, (a, b, order) in enumerate(edges)]
    return {"format": FORMAT, "id": slug, "name": name, "atoms": atoms, "bonds": bonds,
            "ports": [{"atomId": atom["id"], "maxOrder": max_order} for atom in atoms],
            "rigid": rigid, "has3d": True, "fragments": []}


def benzene():
    points = [(1.4 * math.cos(i * math.pi / 3), 1.4 * math.sin(i * math.pi / 3), 0) for i in range(6)]
    return _base("苯环", "benzene", points, [(i, (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)])


def ethynyl():
    return _base("炔基", "ethynyl", [(-0.6, 0, 0), (0.6, 0, 0)], [(0, 1, 3)])


def fullerene():
    phi = (1 + math.sqrt(5)) / 2
    vertices = []
    for first, second in itertools.product((-1, 1), repeat=2):
        vertices.extend(((0, first, second * phi), (first, second * phi, 0), (second * phi, 0, first)))
    verts = np.array(vertices, dtype=float)
    adjacency = {i: [j for j in range(12) if abs(np.linalg.norm(verts[i] - verts[j]) - 2) < 1e-8] for i in range(12)}
    directed = [(i, j) for i in range(12) for j in adjacency[i]]
    lookup = {edge: index for index, edge in enumerate(directed)}
    points = [(2 * verts[i] + verts[j]) / 3 * (1.42 / (2 / 3)) for i, j in directed]
    edges = []
    for index, (i, j) in enumerate(directed):
        reverse = lookup[(j, i)]
        if index < reverse:
            edges.append((index, reverse, 2))
        for k in adjacency[i]:
            second = lookup[(i, k)]
            if index < second and k in adjacency[j]:
                edges.append((index, second, 1))
    # maxOrder is the maximum *allowed attachment bond*, not spare valence.
    # The editor dynamically rejects attachments while actual valence is full.
    template = _base("C60 富勒烯", "c60", points, edges, max_order=1, rigid=True)
    template["description"] = "60 个真实碳原子、90 条键；默认整体移动。每个碳已为价态 4，连接前须明确降低相邻双键键级（或设计真实加成反应），不能直接追加第五价键。"
    template["connectionPolicy"] = "requires-explicit-bond-order-edit"
    # Give the actual spherical conformer a distinct planar graph layout.
    from .service import to_mol
    molecule = to_mol(template, sanitize=True)
    for atom, point in zip(template["atoms"], _layout(molecule)):
        atom["pos2"] = point
    assert len(template["atoms"]) == 60 and len(template["bonds"]) == 90
    assert all(sum(b["order"] for b in template["bonds"] if a["id"] in (b["a"], b["b"])) == 4 for a in template["atoms"])
    return template


def write_assets(root=None):
    root = Path(root) if root else Path(__file__).resolve().parents[2] / "frontend" / "editor"
    root.mkdir(parents=True, exist_ok=True)
    data = {"format": "molecular-car-fragment-library-v1", "source": "Procedurally generated from chemical graphs; see molecular_car/editor/fragments.py",
            "fragments": [benzene(), ethynyl(), fullerene()]}
    (root / "fragments.json").write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    return data


if __name__ == "__main__":
    result = write_assets()
    print([(item["id"], len(item["atoms"]), len(item["bonds"])) for item in result["fragments"]])
