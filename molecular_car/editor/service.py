"""RDKit file, layout, constrained geometry and race-admission helpers.

The caller owns edit revisions and must discard stale asynchronous responses.
This module returns the request revision unchanged. Topology, IDs and fragment
records are never replaced by results of a geometry operation.
"""
from __future__ import annotations

import copy
import math
from pathlib import Path
import re
import uuid

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Geometry import Point3D

from ..molecule import MAX_ATOMS, EMBED_SEED, prepare_molecule
from .superatoms import read_groups, write_groups
from .construction import ConstructionError, construct_new_regions, embed_initial

FORMAT = "molecular-car-workshop-v1"
MAX_EDITOR_ATOMS = 2000
BOND_TYPES = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
              3: Chem.BondType.TRIPLE, 1.5: Chem.BondType.AROMATIC}
CHIRAL_TAGS = {name: getattr(Chem.ChiralType, name) for name in (
    "CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW")}
BOND_STEREO = {name: getattr(Chem.BondStereo, name) for name in (
    "STEREONONE", "STEREOANY", "STEREOZ", "STEREOE", "STEREOCIS", "STEREOTRANS")}
BOND_DIRECTIONS = {name: getattr(Chem.BondDir, name) for name in (
    "NONE", "BEGINWEDGE", "BEGINDASH", "ENDDOWNRIGHT", "ENDUPRIGHT", "EITHERDOUBLE", "UNKNOWN")}


class EditorError(ValueError):
    """Human-readable editor operation error, independent of race state."""


def _uid(prefix):
    return prefix + uuid.uuid4().hex[:16]


def _integer(value, label, low=-100, high=100):
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise EditorError(f"{label} 必须是 {low} 到 {high} 之间的整数。")
    return value


def _coords(value, length, label):
    if not isinstance(value, (list, tuple)) or len(value) != length or any(
        isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or abs(v) > 100000
        for v in value
    ):
        raise EditorError(f"{label} 必须包含 {length} 个有限坐标。")
    return [float(v) for v in value]


def check_project(project):
    """Validate serialization and references, deliberately allowing bad chemistry."""
    if not isinstance(project, dict) or project.get("format") != FORMAT:
        raise EditorError("原生项目格式不正确（需要 molecular-car-workshop-v1）。")
    atoms, bonds = project.get("atoms"), project.get("bonds")
    if not isinstance(atoms, list) or not isinstance(bonds, list):
        raise EditorError("项目需要 atoms 和 bonds 数组。")
    if len(atoms) > MAX_EDITOR_ATOMS or len(bonds) > 6000:
        raise EditorError("此编辑器原型最多处理 2000 个原子和 6000 条键。")
    ids = set()
    for atom in atoms:
        if not isinstance(atom, dict) or not isinstance(atom.get("id"), str) or not atom["id"] or atom["id"] in ids:
            raise EditorError("每个原子必须有唯一且稳定的字符串 ID。")
        ids.add(atom["id"])
        try:
            with rdBase.BlockLogs():
                Chem.Atom(atom.get("element", ""))
        except Exception as exc:
            raise EditorError(f"无法识别元素 {atom.get('element')}。") from exc
        _integer(atom.get("charge", 0), "形式电荷")
        _integer(atom.get("isotope", 0), "同位素质量数", 0, 400)
        _integer(atom.get("explicitHs", 0), "显式氢数", 0, 8)
        _integer(atom.get("radicalElectrons", 0), "未成对电子数", 0, 8)
        _coords(atom.get("pos2"), 2, "二维坐标")
        _coords(atom.get("pos3"), 3, "三维坐标")
        if atom.get("chiralTag", "CHI_UNSPECIFIED") not in CHIRAL_TAGS:
            raise EditorError("当前原型支持四面体 CW/CCW 立体标记。")
    bond_ids, pairs = set(), set()
    for bond in bonds:
        if not isinstance(bond, dict) or not isinstance(bond.get("id"), str) or not bond["id"] or bond["id"] in bond_ids:
            raise EditorError("每条键必须有唯一且稳定的字符串 ID。")
        bond_ids.add(bond["id"])
        a, b = bond.get("a"), bond.get("b")
        if a not in ids or b not in ids or a == b:
            raise EditorError("化学键必须连接两个存在且不同的原子。")
        pair = frozenset((a, b))
        if pair in pairs:
            raise EditorError("两个原子之间不能有重复键；请修改已有键的键级。")
        pairs.add(pair)
        if isinstance(bond.get("order"), bool) or bond.get("order") not in BOND_TYPES:
            raise EditorError("当前文件适配器支持单键、双键、三键和芳香键。")
        if bond.get("stereo", "STEREONONE") not in BOND_STEREO or bond.get("direction", "NONE") not in BOND_DIRECTIONS:
            raise EditorError("存在未支持的键立体标记。")
        stereo_atoms = bond.get("stereoAtoms", [])
        if stereo_atoms and (len(stereo_atoms) != 2 or any(i not in ids for i in stereo_atoms)):
            raise EditorError("立体键的参考原子 ID 无效。请重新指定立体化学。")
    for fragment in project.get("fragments", []):
        if any(atom_id not in ids for atom_id in fragment.get("atomIds", [])):
            raise EditorError("片段引用了不存在的原子。")
    return project


def to_mol(project, dimension=3, sanitize=False):
    check_project(project)
    builder = Chem.RWMol()
    indices = {atom["id"]: index for index, atom in enumerate(project["atoms"])}
    for record in project["atoms"]:
        atom = Chem.Atom(record["element"])
        atom.SetFormalCharge(record.get("charge", 0))
        atom.SetIsotope(record.get("isotope", 0))
        atom.SetNumExplicitHs(record.get("explicitHs", 0))
        atom.SetNoImplicit(bool(record.get("noImplicit", False)))
        atom.SetNumRadicalElectrons(record.get("radicalElectrons", 0))
        atom.SetChiralTag(CHIRAL_TAGS[record.get("chiralTag", "CHI_UNSPECIFIED")])
        builder.AddAtom(atom)
    for record in project["bonds"]:
        builder.AddBond(indices[record["a"]], indices[record["b"]], BOND_TYPES[record["order"]])
        bond = builder.GetBondWithIdx(builder.GetNumBonds() - 1)
        if record["order"] == 1.5:
            bond.SetIsAromatic(True)
            bond.GetBeginAtom().SetIsAromatic(True)
            bond.GetEndAtom().SetIsAromatic(True)
        bond.SetBondDir(BOND_DIRECTIONS[record.get("direction", "NONE")])
    molecule = builder.GetMol()
    molecule.UpdatePropertyCache(strict=False)
    for record, bond in zip(project["bonds"], molecule.GetBonds()):
        refs = record.get("stereoAtoms", [])
        if refs:
            try:
                bond.SetStereoAtoms(indices[refs[0]], indices[refs[1]])
            except Exception as exc:
                raise EditorError("立体键的参考原子不再与该键相邻；请重新指定 E/Z。") from exc
        bond.SetStereo(BOND_STEREO[record.get("stereo", "STEREONONE")])
    conformer = Chem.Conformer(len(project["atoms"]))
    conformer.Set3D(dimension == 3)
    for index, atom in enumerate(project["atoms"]):
        position = atom["pos3"] if dimension == 3 else [*atom["pos2"], 0]
        conformer.SetAtomPosition(index, Point3D(*position))
    molecule.AddConformer(conformer)
    molecule.SetProp("_Name", str(project.get("name", "Molecular car"))[:120].replace("\n", " ").replace("\r", " "))
    if sanitize:
        try:
            with rdBase.BlockLogs():
                Chem.SanitizeMol(molecule)
        except Exception as exc:
            raise EditorError(f"无法计算几何。请检查价态、形式电荷和芳香环。({exc})") from exc
    return molecule


def _depiction_copy(molecule):
    """Perceive drawing geometry without sanitizing or changing the graph.

    Layout accepts unfinished drawings but needs hybridization to depict sp
    carbons at 180 degrees. Perception happens on an indexed private copy;
    canonical bond orders and the player's 3D conformer remain untouched.
    """
    layout = Chem.Mol(molecule)
    layout.UpdatePropertyCache(strict=False)
    Chem.SetHybridization(layout)
    return layout


def _layout(molecule):
    try:
        with rdBase.BlockLogs():
            layout = _depiction_copy(molecule)
            rdDepictor.Compute2DCoords(layout, clearConfs=True)
        return layout.GetConformer().GetPositions()[:, :2].tolist()
    except Exception:
        # Incomplete graphs are still editable, even if depiction cannot lay out.
        return [[index * 1.5, 0.0] for index in range(molecule.GetNumAtoms())]


def _resolve_drawn_stereo(project):
    """Resolve user-drawn 2D wedge/dash declarations before using any conformer.

    A wedge is a declaration relative to the 2D layout, including when has3d is
    already true. Merely SetBondDir is insufficient for ETKDG. Work on indexed
    copies, reject contradictions, and persist newly resolved tags in the output
    snapshot of an explicit geometry operation.
    """
    drawn = [b for b in project["bonds"] if b.get("direction") in ("BEGINWEDGE", "BEGINDASH")]
    if not drawn:
        return []
    two_d = to_mol(project, 2)
    for atom in two_d.GetAtoms():
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    with rdBase.BlockLogs():
        Chem.AssignChiralTypesFromBondDirs(two_d, replaceExistingTags=True)
    index_by_id = {a["id"]: i for i, a in enumerate(project["atoms"])}
    centers = {index_by_id[b["a"]] for b in drawn}
    warnings = []
    for index in centers:
        inferred = str(two_d.GetAtomWithIdx(index).GetChiralTag())
        declared = project["atoms"][index].get("chiralTag", "CHI_UNSPECIFIED")
        if inferred == "CHI_UNSPECIFIED":
            raise EditorError(f"原子 {index + 1} 的二维楔线无法确定四面体手性；请从立体中心指向邻接原子画楔线，并避免共线/重叠布局。")
        if declared != "CHI_UNSPECIFIED" and inferred != declared:
            raise EditorError(f"原子 {index + 1} 的二维楔/虚楔与已保存的手性标记冲突。请修正二维楔线或更新手性标记后重试。")
        if declared == "CHI_UNSPECIFIED":
            project["atoms"][index]["chiralTag"] = inferred
            warnings.append(f"已识别原子 {index + 1} 的手性。")
    return warnings


def _apply_stereo_layout(project, molecule):
    """Change 2D depiction and regenerate wedge orientation for the same tags."""
    with rdBase.BlockLogs():
        laid_out = _depiction_copy(molecule)
        rdDepictor.Compute2DCoords(laid_out, clearConfs=True)
        for bond in laid_out.GetBonds():
            if bond.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH):
                bond.SetBondDir(Chem.BondDir.NONE)
        Chem.WedgeMolBonds(laid_out, laid_out.GetConformer())
    for atom, position in zip(project["atoms"], laid_out.GetConformer().GetPositions()):
        atom["pos2"] = position[:2].tolist()
    for record, bond in zip(project["bonds"], laid_out.GetBonds()):
        record["direction"] = str(bond.GetBondDir())
        # RDKit can reverse the stored endpoint order to draw from the center.
        # Graph connections and bond indices remain unchanged.
        record["a"] = project["atoms"][bond.GetBeginAtomIdx()]["id"]
        record["b"] = project["atoms"][bond.GetEndAtomIdx()]["id"]


def _project_from_mol(molecule, filename):
    has_3d = bool(molecule.GetNumConformers() and molecule.GetConformer().Is3D())
    points = molecule.GetConformer().GetPositions() if molecule.GetNumConformers() else np.zeros((molecule.GetNumAtoms(), 3))
    pos2 = _layout(molecule) if has_3d else points[:, :2].tolist()
    atom_ids = [_uid("a-") for _ in molecule.GetAtoms()]
    atoms = []
    for index, atom in enumerate(molecule.GetAtoms()):
        atoms.append({"id": atom_ids[index], "element": atom.GetSymbol(), "charge": atom.GetFormalCharge(),
                      "isotope": atom.GetIsotope(), "explicitHs": atom.GetNumExplicitHs(),
                      "noImplicit": atom.GetNoImplicit(), "radicalElectrons": atom.GetNumRadicalElectrons(),
                      "chiralTag": str(atom.GetChiralTag()), "pos2": pos2[index], "pos3": points[index].tolist()})
    bonds = [{"id": _uid("b-"), "a": atom_ids[bond.GetBeginAtomIdx()], "b": atom_ids[bond.GetEndAtomIdx()],
              "order": bond.GetBondTypeAsDouble(), "stereo": str(bond.GetStereo()),
              "stereoAtoms": [atom_ids[i] for i in bond.GetStereoAtoms()], "direction": str(bond.GetBondDir())}
             for bond in molecule.GetBonds()]
    project = {"format": FORMAT, "name": molecule.GetProp("_Name") if molecule.HasProp("_Name") else Path(filename).stem,
               "atoms": atoms, "bonds": bonds, "fragments": [], "has3d": has_3d,
               "metadata": {"sourceFilename": filename, "rdkitVersion": rdBase.rdkitVersion,
                            "indexToAtomId": atom_ids}}
    groups = read_groups(molecule, atom_ids)
    if groups:
        project["metadata"]["ketcherGroups"] = groups
    if has_3d:
        # A 3D MOL parser assigns chirality from coordinates and normally drops
        # its drawing wedges. The new 2D depiction needs its own correctly
        # oriented wedges; a raw CW/CCW tag is not a Ketcher stereoParity value.
        # This works on a molecule copy and transfers only indexed depiction
        # fields. Input pos3, stable IDs and chiral tags remain untouched.
        try:
            depicted = copy.deepcopy(project)
            _apply_stereo_layout(depicted, molecule)
            project = depicted
        except Exception as exc:
            # Incomplete molecules remain importable and saveable even when
            # RDKit cannot yet produce a chemically consistent depiction.
            project["metadata"]["stereoNotice"] = (
                "未能生成完整的二维立体描画；"
                f"请修正结构后整理二维排版。({exc})"
            )
    return project


def import_file(data):
    filename, content = data.get("filename", "molecule.mol"), data.get("content")
    if not isinstance(filename, str) or Path(filename).suffix.lower() not in {".mol", ".sdf"}:
        raise EditorError("请导入 MOL 或单分子 SDF 文件。原生项目由编辑器直接打开。")
    if not isinstance(content, str) or not content.strip() or len(content) > 2_000_000:
        raise EditorError("分子文件为空或超过 2 MB。")
    normalized = content.removeprefix("\ufeff").replace("\r\n", "\n")
    records = [record for record in re.split(r"^\$\$\$\$[ \t]*$", normalized, flags=re.MULTILINE) if record.strip()]
    if len(records) != 1:
        raise EditorError(f"此 SDF 包含 {len(records)} 个分子记录。当前原型只打开单记录文件；请先选择一个记录单独导出。")
    if len(re.findall(r"^M\s+END\s*$", records[0], flags=re.MULTILINE)) != 1:
        raise EditorError("MOL/SDF 的结构记录不完整，或包含多个结构。")
    try:
        with rdBase.BlockLogs():
            molecule = Chem.MolFromMolBlock(records[0], sanitize=False, removeHs=False, strictParsing=True)
        if molecule is None:
            raise EditorError("无法读取原子、键或坐标记录。")
        if molecule.GetNumAtoms() > MAX_EDITOR_ATOMS:
            raise EditorError("此编辑器原型最多处理 2000 个原子。")
        if any(atom.HasQuery() for atom in molecule.GetAtoms()) or any(bond.HasQuery() for bond in molecule.GetBonds()):
            raise EditorError("查询原子和查询键尚未支持；请导出具体元素与常规键。")
        if molecule.GetStereoGroups() or any(group.GetProp("TYPE") != "SUP" for group in Chem.GetMolSubstanceGroups(molecule)):
            raise EditorError("此文件包含增强立体组或非 Superatom 的 S-group；当前编辑器暂不支持此类型。")
        warnings = []
        cleaned = Chem.Mol(molecule)
        try:
            with rdBase.BlockLogs():
                Chem.SanitizeMol(cleaned)
                # A 3D MOL parser can attach a geometrical CW/CCW tag even to
                # CH3/NH3 centers with equivalent substituents. These are not
                # stereocenters; canonical stereo perception removes them.
                Chem.AssignStereochemistry(cleaned, cleanIt=True, force=True)
            # Keep the input bond orders and atom order; only transfer stereo
            # assigned from input coordinates/wedges on this exact indexed graph.
            for raw, clean in zip(molecule.GetAtoms(), cleaned.GetAtoms()):
                raw.SetChiralTag(clean.GetChiralTag())
            for raw, clean in zip(molecule.GetBonds(), cleaned.GetBonds()):
                if len(clean.GetStereoAtoms()) == 2:
                    raw.SetStereoAtoms(*clean.GetStereoAtoms())
                raw.SetStereo(clean.GetStereo())
        except Exception as exc:
            warnings.append(f"结构未通过化学检查：{exc}")
        project = _project_from_mol(molecule, filename)
        check_project(project)
        if project["metadata"].get("stereoNotice"):
            warnings.append(project["metadata"]["stereoNotice"])
        return {"project": project, "warnings": warnings,
                "mapping": {"indexToAtomId": [a["id"] for a in project["atoms"]]}}
    except EditorError:
        raise
    except Exception as exc:
        raise EditorError(f"读取分子文件失败：{exc}") from exc


def _stereo_signature(molecule):
    copy_mol = Chem.Mol(molecule)
    with rdBase.BlockLogs():
        Chem.AssignStereochemistry(copy_mol, cleanIt=True, force=True)
    atoms = {atom.GetIdx(): atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else str(atom.GetChiralTag())
             for atom in copy_mol.GetAtoms() if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED}
    bonds = {bond.GetIdx(): str(bond.GetStereo()).replace("STEREOTRANS", "STEREOE").replace("STEREOCIS", "STEREOZ")
             for bond in copy_mol.GetBonds() if bond.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)}
    return atoms, bonds


def _check_stereo_roundtrip(molecule, block, dimension):
    """Index-based stereo check: no graph correspondence or SMILES identity."""
    expected_atoms, expected_bonds = _stereo_signature(molecule)
    with rdBase.BlockLogs():
        observed = Chem.MolFromMolBlock(block, sanitize=True, removeHs=False)
    if observed is None:
        raise EditorError("导出文件校验失败。请检查结构后重试。")
    got_atoms, got_bonds = _stereo_signature(observed)
    errors = []
    for atom in molecule.GetAtoms():
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED and atom.GetIdx() not in expected_atoms:
            errors.append(f"原子 {atom.GetIdx() + 1} 的手性标记不对应有效四面体立体中心")
    for bond in molecule.GetBonds():
        if bond.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY) and bond.GetIdx() not in expected_bonds:
            errors.append(f"键 {bond.GetIdx() + 1} 的立体标记不对应有效双键立体中心")
    for index, signature in expected_atoms.items():
        if got_atoms.get(index) != signature:
            errors.append(f"原子 {index + 1} 的手性标记与坐标不一致")
    for index, signature in expected_bonds.items():
        if got_bonds.get(index) != signature:
            errors.append(f"键 {index + 1} 的 E/Z 标记与坐标不一致")
    if errors:
        action = "生成符合标记的三维构象，或有意修改立体标记" if dimension == 3 else "整理二维排版以满足 E/Z 标记"
        raise EditorError("；".join(errors) + f"。请先{action}。")
    warnings = []
    if dimension == 3 and (set(got_atoms) - set(expected_atoms) or set(got_bonds) - set(expected_bonds)):
        warnings.append("3D MOL 会使外部软件从坐标推定当前尚未显式标记的立体化学。")
    return warnings


def export_file(data):
    project = copy.deepcopy(check_project(data.get("project")))
    dimension, file_format = data.get("dimension", 3), data.get("format", "mol")
    if dimension not in (2, 3) or file_format not in ("mol", "sdf"):
        raise EditorError("导出维度应为 2 或 3，文件格式应为 mol 或 sdf。")
    if dimension == 3 and not project.get("has3d"):
        raise EditorError("此项目尚未生成三维构象。请先生成三维，或导出二维 MOL。")
    warnings = _resolve_drawn_stereo(project)
    molecule = to_mol(project, dimension)
    valid_chemistry = True
    try:
        with rdBase.BlockLogs():
            Chem.SanitizeMol(Chem.Mol(molecule))
    except Exception as exc:
        valid_chemistry = False
        warnings.append(f"结构未通过化学检查，其他软件可能无法读取此文件。建议同时保存工程。({exc})")
    try:
        write_groups(molecule, project)
        with rdBase.BlockLogs():
            # Do not kekulize: preserve explicit single/double/aromatic orders.
            block = Chem.MolToMolBlock(molecule, includeStereo=True, kekulize=False,
                                      forceV3000=len(project["atoms"]) > 999 or len(project["bonds"]) > 999
                                      or bool(project.get("metadata", {}).get("ketcherGroups")))
        if valid_chemistry:
            warnings.extend(_check_stereo_roundtrip(molecule, block, dimension))
        elif any(a.get("chiralTag", "CHI_UNSPECIFIED") != "CHI_UNSPECIFIED" for a in project["atoms"]):
            warnings.append("当前价态未通过检查，无法认证导出后的立体信息。建议同时保存工程。")
    except EditorError:
        raise
    except Exception as exc:
        raise EditorError(f"无法导出此结构。({exc})") from exc
    if file_format == "sdf":
        block += "$$$$\n"
    name = re.sub(r"[<>:\"/\\|?*\x00-\x1f]", "_", str(project.get("name", "my-car")))[:80] or "my-car"
    return {"content": block, "filename": f"{name}.{file_format}", "warnings": warnings,
            "mapping": {"indexToAtomId": [a["id"] for a in project["atoms"]]}}


def _materialize_hydrogens(project, work, original_count):
    """Append AddHs results by index without reassigning any existing identity.

    All old pos2/pos3 values are untouched here. Generate's caller explicitly
    updates old pos3 before invoking this helper; hydrogenate never does so.
    """
    added_count = work.GetNumAtoms() - original_count
    if not added_count:
        return {"hydrogensAdded": 0, "addedAtomIds": [], "addedBondIds": []}
    if work.GetNumAtoms() > MAX_EDITOR_ATOMS:
        raise EditorError(f"补氢后将有 {work.GetNumAtoms()} 个原子，超过编辑器的 {MAX_EDITOR_ATOMS} 个原子上限。")
    # AddHs places only newly appended H around the supplied 2D conformer.
    # This preserves the player's old layout without a second graph matching
    # or depiction pass. Both AddHs calls start from the same indexed graph.
    with rdBase.BlockLogs():
        drawing = Chem.AddHs(to_mol(project, 2, sanitize=True), addCoords=True)
    if drawing.GetNumAtoms() != work.GetNumAtoms() or any(
        drawing.GetAtomWithIdx(i).GetAtomicNum() != work.GetAtomWithIdx(i).GetAtomicNum()
        for i in range(work.GetNumAtoms())
    ):
        raise EditorError("补氢的二维/三维编号映射不一致。")
    points2 = drawing.GetConformer().GetPositions()
    points3 = work.GetConformer().GetPositions()
    if not np.isfinite(points2).all() or not np.isfinite(points3).all():
        raise EditorError("当前局部坐标无法放置氢原子；请先修正重叠原子或生成三维。")
    original_bond_count = len(project["bonds"])
    ids = [atom["id"] for atom in project["atoms"]]
    added_ids = [_uid("a-") for _ in range(added_count)]
    ids.extend(added_ids)
    parents = {}
    for index in range(original_count, work.GetNumAtoms()):
        atom = work.GetAtomWithIdx(index)
        if atom.GetAtomicNum() != 1 or atom.GetDegree() != 1:
            raise EditorError("补氢返回了非预期原子或连接。")
        parent = atom.GetNeighbors()[0].GetIdx()
        parents[ids[index]] = ids[parent]
        project["atoms"].append({"id": ids[index], "element": "H", "charge": atom.GetFormalCharge(),
                                 "isotope": atom.GetIsotope(), "explicitHs": 0,
                                 "noImplicit": atom.GetNoImplicit(), "radicalElectrons": atom.GetNumRadicalElectrons(),
                                 "chiralTag": str(atom.GetChiralTag()),
                                 "pos2": points2[index, :2].tolist(), "pos3": points3[index].tolist(),
                                 "generatedHydrogen": True})
    for index, record in enumerate(project["atoms"][:original_count]):
        atom = work.GetAtomWithIdx(index)
        # Bracket H counts (e.g. [C@H], [NH3+]) have now become real graph
        # vertices. Retaining the old count would double-count those hydrogens.
        if record.get("explicitHs", 0) != atom.GetNumExplicitHs():
            record["explicitHs"] = atom.GetNumExplicitHs()
        if record.get("chiralTag", "CHI_UNSPECIFIED") != str(atom.GetChiralTag()):
            record["chiralTag"] = str(atom.GetChiralTag())
    added_bond_ids = []
    for bond in list(work.GetBonds())[original_bond_count:]:
        new_id = _uid("b-")
        added_bond_ids.append(new_id)
        project["bonds"].append({"id": new_id, "a": ids[bond.GetBeginAtomIdx()], "b": ids[bond.GetEndAtomIdx()],
                                 "order": bond.GetBondTypeAsDouble(), "stereo": str(bond.GetStereo()),
                                 "stereoAtoms": [ids[i] for i in bond.GetStereoAtoms()], "direction": "NONE"})
    # New H belongs to every fragment containing its parent; rigid fragments
    # must carry these actual atoms on the next player drag as well.
    for fragment in project.get("fragments", []):
        members = set(fragment["atomIds"])
        fragment["atomIds"].extend(atom_id for atom_id in added_ids if parents[atom_id] in members)
    # Replacing an implicit H by a plotted fourth neighbor changes how a 2D
    # wedge is interpreted. Draw the newly materialized H consistently while
    # keeping all pre-existing coordinates, bond IDs and absolute stereo.
    centers = {parents[atom_id] for atom_id in added_ids}
    centers = {i for i in range(original_count) if ids[i] in centers
               and work.GetAtomWithIdx(i).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED}
    if centers:
        depiction = to_mol(project, 2, sanitize=True)
        for bond in depiction.GetBonds():
            if bond.GetBeginAtomIdx() in centers and bond.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH):
                bond.SetBondDir(Chem.BondDir.NONE)
        for index in centers:
            new_h = next(a.GetIdx() for a in depiction.GetAtomWithIdx(index).GetNeighbors()
                         if a.GetIdx() >= original_count and a.GetAtomicNum() == 1)
            Chem.WedgeBond(depiction.GetBondBetweenAtoms(index, new_h), index, depiction.GetConformer())
        for record, bond in zip(project["bonds"], depiction.GetBonds()):
            if record.get("direction", "NONE") != str(bond.GetBondDir()):
                record["direction"] = str(bond.GetBondDir())
            record["a"], record["b"] = ids[bond.GetBeginAtomIdx()], ids[bond.GetEndAtomIdx()]
    check_project(project)
    return {"hydrogensAdded": added_count, "addedAtomIds": added_ids, "addedBondIds": added_bond_ids}


def _pending_geometry(project):
    """Explicit edit provenance, never inferred from distances or graph matching."""
    pending = project.get("metadata", {}).get("geometryPending")
    if pending is None:
        return set(), set(), set()
    if not isinstance(pending, dict):
        raise EditorError("待准备的三维编辑记录必须是对象。")
    atom_ids = {atom["id"] for atom in project["atoms"]}
    result = []
    for name in ("atomIds", "changedAtomIds", "editedAtomIds"):
        values = pending.get(name, [])
        if not isinstance(values, list) or any(not isinstance(value, str) or value not in atom_ids for value in values):
            raise EditorError("待准备的三维编辑记录包含不存在的原子 ID；请重新编辑此区域。")
        result.append(set(values))
    return tuple(result)


def _c60_geometry_seeds(project, new_ids, first_generation):
    """Use the native template's explicit index mapping and genuine sphere seed."""
    groups = project.get("metadata", {}).get("ketcherGroups", [])
    atoms = {atom["id"]: atom for atom in project["atoms"]}
    seeds = []
    for group in groups if isinstance(groups, list) else []:
        if not isinstance(group, dict) or group.get("kind") != "c60":
            continue
        members = group.get("atomIds", [])
        if not isinstance(members, list) or any(not isinstance(member, str) for member in members):
            raise EditorError("C60 原生模板的稳定原子映射不完整；请展开并检查结构。")
        if not first_generation and not new_ids.intersection(members):
            continue
        if (not isinstance(members, list) or len(members) != 60 or len(set(members)) != 60
                or any(member not in atoms or atoms[member]["element"] != "C" for member in members)):
            raise EditorError("C60 原生模板的稳定原子映射不完整；请展开并检查结构。")
        points = np.array([atoms[member]["pos3"] for member in members], dtype=float)
        centered = points - points.mean(axis=0)
        singular = np.linalg.svd(centered, compute_uv=False)
        internal = [bond for bond in project["bonds"] if bond["a"] in members and bond["b"] in members]
        lengths = [np.linalg.norm(np.array(atoms[bond["a"]]["pos3"]) - atoms[bond["b"]]["pos3"])
                   for bond in internal]
        if (len(internal) != 90 or singular[-1] < 0.1 or not lengths
                or min(lengths) < 0.8 or max(lengths) > 2.0):
            raise EditorError("C₆₀ 球壳坐标无效。请重新添加 C₆₀ 模板。")
        seeds.append((members, points))
    return seeds


def _fixed_pending_connection_warnings(project, changed_ids, fixed_ids):
    """Warn about this edit's strained fixed boundary, never scan unchanged cars.

    Covalent radii provide a deliberately broad advisory range rather than a
    chemistry acceptance rule. Automatic construction leaves both existing
    endpoints intact even if their input connection length is unrealistic.
    """
    if not changed_ids:
        return []
    atoms = {atom["id"]: atom for atom in project["atoms"]}
    periodic_table = Chem.GetPeriodicTable()
    abnormal = []
    for bond in project["bonds"]:
        endpoints = {bond["a"], bond["b"]}
        if not endpoints.issubset(fixed_ids) or not endpoints.intersection(changed_ids):
            continue
        a, b = atoms[bond["a"]], atoms[bond["b"]]
        if a["element"] == "H" or b["element"] == "H":
            continue
        radius_sum = sum(periodic_table.GetRcovalent(atom["element"]) for atom in (a, b))
        length = float(np.linalg.norm(np.array(a["pos3"]) - b["pos3"]))
        if radius_sum > 0 and not 0.5 * radius_sum <= length <= 1.8 * radius_sum:
            abnormal.append(f"{a['element']}–{b['element']} {length:.2f} Å")
    if not abnormal:
        return []
    examples = "、".join(abnormal[:3]) + ("等" if len(abnormal) > 3 else "")
    return [f"连接键长异常（{examples}）。请移动片段或运行 UFF 优化。"]


def _prepare_geometry(data, project, *, rebuild=False):
    """Initial coordinates and real H only; automatic switching never uses UFF."""
    first_generation = rebuild or not project.get("has3d")
    new_ids, changed_ids, edited_ids = _pending_geometry(project)
    seeds = _c60_geometry_seeds(project, new_ids, first_generation)
    if not first_generation and not (new_ids or changed_ids or edited_ids):
        result = geometry({**data, "project": project, "operation": "hydrogenate"})
        result["construction"] = {"optimized": False, "method": "hydrogenate-only",
                                  "existingAtomsPreserved": len(project["atoms"]), "regions": []}
    else:
        warnings = _resolve_drawn_stereo(project)
        molecule = to_mol(project, 2 if first_generation and not seeds else 3, sanitize=True)
        count = molecule.GetNumAtoms()
        if count > 600:
            raise EditorError("三维构建最多处理 600 个原子。请减少原子数量后重试。")
        ids = {atom["id"]: index for index, atom in enumerate(project["atoms"])}
        fixed_ids = data.get("fixedIds", [])
        if not isinstance(fixed_ids, list) or any(value not in ids for value in fixed_ids):
            raise EditorError("固定原子列表包含无效 ID。")
        if first_generation and fixed_ids:
            raise EditorError("固定原子需要已有三维坐标；请先生成三维。")
        # Changes to old chemistry do not grant permission to move old atoms.
        # Their exact coordinates, including every existing H, remain intact.
        movable = set(range(count)) if first_generation else {ids[value] for value in new_ids}
        explicit_fixed = {ids[value] for value in fixed_ids}
        movable.difference_update(explicit_fixed)
        fixed = set(range(count)) - movable
        try:
            if first_generation and not seeds:
                with rdBase.BlockLogs():
                    work = Chem.AddHs(molecule)
                embed_initial(work, EMBED_SEED)
                regions = [{"atomCount": count, "method": "initial-ETKDG"}]
            else:
                seed_groups = [{ids[value] for value in members} for members, _ in seeds]
                cages = [{ids[value] for value in group.get("atomIds", []) if value in ids}
                         for group in project.get("metadata", {}).get("ketcherGroups", [])
                         if isinstance(group, dict) and group.get("kind") == "c60"]
                regions = construct_new_regions(molecule, movable, seed_groups, cages, explicit_fixed, EMBED_SEED)
                # AddHs sees the finished initial heavy-atom positions, never
                # provisional 2D units or an embedded replacement of the car.
                with rdBase.BlockLogs():
                    work = Chem.AddHs(molecule, addCoords=True)
        except ConstructionError as exc:
            raise EditorError(str(exc)) from exc
        if ((work.HasProp("_workshopCurvedFallback") and work.GetBoolProp("_workshopCurvedFallback"))
                or any(region.get("curvedFallback") for region in regions)):
            warnings.append("已生成弯曲构象，可用 UFF 优化改善局部几何。")
        points = work.GetConformer().GetPositions()
        if not np.isfinite(points).all():
            raise EditorError("初始构象含无效坐标。")
        for index in movable:
            project["atoms"][index]["pos3"] = points[index].tolist()
        project["has3d"] = True
        # Rebuilding coordinates uses H in a private construction copy only.
        # The three explicit UI operations must not silently change each other.
        hydrogen_info = (_materialize_hydrogens(project, work, count) if not rebuild
                         else {"hydrogensAdded": 0, "addedAtomIds": [], "addedBondIds": []})
        warnings.extend(export_file({"project": project, "dimension": 3, "format": "mol"})["warnings"])
        warnings.extend(_fixed_pending_connection_warnings(
            project, changed_ids, {project["atoms"][index]["id"] for index in fixed}))
        message = ("三维构象已重建。" if rebuild else "三维结构已就绪。")
        if hydrogen_info["hydrogensAdded"]:
            message = message.rstrip("。") + f" · 已添加 {hydrogen_info['hydrogensAdded']} 个氢原子。"
        result = {"project": project, "revision": data.get("revision"), "message": message, "warnings": warnings, **hydrogen_info,
                  "construction": {"optimized": False, "method": "initial-ETKDG" if first_generation and not seeds else "local-construction",
                                   "existingAtomsPreserved": len(fixed), "regions": regions},
                  "mapping": {"indexToAtomId": [atom["id"] for atom in project["atoms"]]}}
    result["project"].get("metadata", {}).pop("geometryPending", None)
    result["project"].get("metadata", {}).pop("geometryNotice", None)
    return result


def geometry(data):
    project = copy.deepcopy(check_project(data.get("project")))
    operation = data.get("operation")
    if operation not in ("generate", "rebuild", "hydrogenate", "relax", "layout2d", "prepare"):
        raise EditorError("未知的几何操作。")
    if not project["atoms"]:
        raise EditorError("请先添加原子。")
    if operation in ("prepare", "rebuild"):
        return _prepare_geometry(data, project, rebuild=operation == "rebuild")
    warnings = _resolve_drawn_stereo(project)
    molecule = to_mol(project, 3 if project.get("has3d") else 2, sanitize=operation != "layout2d")
    if operation == "layout2d":
        _apply_stereo_layout(project, molecule)
        message = "二维排版已整理。"
    else:
        fixed_ids = data.get("fixedIds", [])
        if not isinstance(fixed_ids, list) or any(i not in {a["id"] for a in project["atoms"]} for i in fixed_ids):
            raise EditorError("固定原子列表包含无效 ID。")
        if operation == "relax" and not project.get("has3d"):
            raise EditorError("请先生成三维构象，再进行几何整理或仅补氢。")
        if fixed_ids and not project.get("has3d"):
            raise EditorError("固定原子需要已有三维坐标；请先生成三维。")
        original_count = molecule.GetNumAtoms()
        if original_count > 600:
            raise EditorError("几何操作最多处理 600 个原子。请减少原子数量后重试。")
        with rdBase.BlockLogs():
            work = Chem.AddHs(molecule, addCoords=operation in ("relax", "hydrogenate"))
        if operation == "hydrogenate":
            hydrogen_info = _materialize_hydrogens(project, work, original_count)
            warnings.extend(export_file({"project": project, "dimension": 3 if project.get("has3d") else 2, "format": "mol"})["warnings"])
            message = (f"已添加 {hydrogen_info['hydrogensAdded']} 个氢原子。" if hydrogen_info["hydrogensAdded"] else "氢原子已补齐。")
            return {"project": project, "revision": data.get("revision"), "message": message, "warnings": warnings, **hydrogen_info,
                    "mapping": {"indexToAtomId": [a["id"] for a in project["atoms"]]}}
        with rdBase.BlockLogs():
            has_parameters = AllChem.UFFHasAllMoleculeParams(work)
        if not has_parameters:
            raise EditorError("此结构缺少 UFF 参数。请手动调整构象。")
        fixed = [index for index, atom in enumerate(project["atoms"]) if atom["id"] in fixed_ids]
        if operation == "generate":
            params = AllChem.ETKDGv3()
            params.randomSeed, params.numThreads = EMBED_SEED, 1
            params.maxIterations = 200
            params.timeout = 20
            if fixed:
                params.SetCoordMap({index: Point3D(*project["atoms"][index]["pos3"]) for index in fixed})
                params.useRandomCoords = True
            with rdBase.BlockLogs():
                status = AllChem.EmbedMolecule(work, params)
            if status != 0 and not fixed:
                params.useRandomCoords = True
                with rdBase.BlockLogs():
                    status = AllChem.EmbedMolecule(work, params)
            if status != 0:
                raise EditorError("无法生成符合结构或固定点约束的三维构象。请检查结构及固定原子。")
            for index in fixed:
                work.GetConformer().SetAtomPosition(index, Point3D(*project["atoms"][index]["pos3"]))
        with rdBase.BlockLogs():
            forcefield = AllChem.UFFGetMoleculeForceField(work, vdwThresh=10.0, ignoreInterfragInteractions=True)
            for index in fixed:
                forcefield.AddFixedPoint(index)
            forcefield.Initialize()
            initial_energy = float(forcefield.CalcEnergy())
            result = forcefield.Minimize(maxIts=1000, forceTol=0.001, energyTol=1e-7)
            final_energy = float(forcefield.CalcEnergy())
        points = work.GetConformer().GetPositions()[:original_count]
        if not np.isfinite(points).all() or not math.isfinite(final_energy):
            raise EditorError("优化产生无效坐标。请重建三维构象后重试。")
        for atom, point in zip(project["atoms"], points):
            atom["pos3"] = point.tolist()
        project["has3d"] = True
        hydrogen_info = (_materialize_hydrogens(project, work, original_count) if operation == "generate"
                         else {"hydrogensAdded": 0, "addedAtomIds": [], "addedBondIds": []})
        # Check geometrical inversion before returning. Original stereo fields
        # remain untouched; an incompatible result is never applied silently.
        warnings.extend(export_file({"project": project, "dimension": 3, "format": "mol"})["warnings"])
        if result:
            warnings.append("UFF 优化尚未收敛（1000 步）。请继续优化或调整构象。")
        if operation == "generate":
            message = f"三维构象已生成 · 已添加 {hydrogen_info['hydrogensAdded']} 个氢原子。"
        else:
            message = "UFF 优化完成。"
        if fixed:
            message = message.rstrip("。") + f" · 固定 {len(fixed)} 个原子。"
        if len(Chem.GetMolFrags(molecule)) > 1:
            warnings.append("当前有多个未连接片段，优化按片段分别进行。请连接后再检查整体构象。")
        return {"project": project, "revision": data.get("revision"), "message": message, "warnings": warnings, **hydrogen_info,
                "energy": {"before": initial_energy, "after": final_energy, "unit": "kcal/mol", "converged": result == 0},
                "mapping": {"indexToAtomId": [a["id"] for a in project["atoms"]]}}
    return {"project": project, "revision": data.get("revision"), "message": message, "warnings": warnings,
            "mapping": {"indexToAtomId": [a["id"] for a in project["atoms"]]}}


def validate(data):
    project = check_project(data.get("project"))
    response = {"editable": True, "eligible": False, "errors": [], "warnings": [],
                "atomCount": len(project["atoms"]), "maxRaceAtoms": MAX_ATOMS,
                "revision": data.get("revision"),
                "preparationChanges": "赛道后端会补氢、UFF 松弛并平移，比赛构象可能与编辑器不同。"}
    try:
        # Reuse exactly the current race preparation on a private molecule. It
        # neither accesses GameController nor starts or contacts LAMMPS.
        exported = export_file({"project": project, "dimension": 3 if project.get("has3d") else 2, "format": "mol"})
        response["warnings"].extend(exported["warnings"])
        prepared = prepare_molecule(exported["content"], exported["filename"])
        span = np.ptp(prepared.positions, axis=0)
        # engine.py performs this final placement constraint, independently of
        # the force-field preparation module. Keep the editor check aligned.
        if float(span.max()) > 45:
            raise EditorError("分子的空间跨度超过赛道的 45 Å 限制。请调整构象或使用尺寸更小的分子。")
        response.update(eligible=True, preparedAtomCount=len(prepared.atoms), span=span.tolist())
        response["warnings"].extend(prepared.notes)
    except (ValueError, RuntimeError) as exc:
        response["errors"].append(str(exc))
    return response


def dispatch(path, data):
    """Return JSON-safe data; raise ValueError for an HTTP 400 response.

    Integration: before /api/load, if path.startswith('/api/editor/'):
        self._json(dispatch(path, data)); return
    Routes do not change or inspect game state. No service lifecycle calls.
    """
    handlers = {"/api/editor/import": import_file, "/api/editor/export": export_file,
                "/api/editor/geometry": geometry, "/api/editor/validate": validate}
    if path not in handlers:
        raise EditorError("未知的编辑器接口。")
    if not isinstance(data, dict):
        raise EditorError("请求必须为 JSON 对象。")
    try:
        return handlers[path](data)
    except EditorError:
        raise
    except Exception as exc:
        # All adapter and native RDKit failures become the host's JSON 400
        # response, including malformed nested fields. Never break HTTP midway.
        raise EditorError(f"编辑器操作未完成：{exc}") from exc
