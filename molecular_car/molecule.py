"""Prepare a chemically explicit molecule for the nanocar prototype.

The RDKit UFF object provides an energy (kcal/mol) and Cartesian gradient
(kcal/mol/angstrom) for LAMMPS's external-force callback.  UFF's bonded,
angle, torsion, inversion and nonbonded terms remain in RDKit; no ad-hoc
spring network replaces the molecular force field.  RDKit UFF itself does
not add Gasteiger electrostatic interactions.  Those charges are used by
the prototype's separately defined external-field model.

Polarizabilities below are fixed *effective game-model coefficients* in
arbitrary model units.  They are approximate atomwise scales, not a validated
molecular polarizability or a predictor of experimental STM response.
They do not change during a race.  UFF and Gasteiger applicability are
checked independently; a missing parameter never becomes a zero force.

Example files are generated from these explicit molecular graphs using
RDKit ETKDGv3 (seed 20261011) followed by UFF minimization:
  benzene: c1ccccc1 (C6H6, 12 atoms with explicit hydrogen)
  1,2,4,5-tetraphenylbenzene:
    c1(-c2ccccc2)c(-c2ccccc2)cc(-c2ccccc2)c(-c2ccccc2)c1
    (C30H22, 52 atoms with explicit hydrogen).
The latter is a four-arm *test molecule*, not an experimental nanocar.
The files are generated locally from the listed graphs, not downloaded
experimental structures.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Any

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from .element_models import EFFECTIVE_POLARIZABILITY
from .start_pose import PRESERVE_START_POSE, place_start_pose


MAX_ATOMS = 256
EMBED_SEED = 20261011
STARTING_MIN_Z = 2.9

class MoleculePreparationError(ValueError):
    """A user molecule cannot be prepared under this model's rules."""


@dataclass
class PreparedMolecule:
    name: str
    atoms: list[dict]
    bonds: list[dict]
    positions: np.ndarray
    forcefield: Any
    notes: list[str]
    # Keep the conformer alive: the UFF object stores pointers to its positions.
    _rdkit_molecule: Any = field(repr=False)
    forcefield_validation: dict = field(default_factory=dict)
    start_pose: dict = field(default_factory=dict)


def _check_uff_cartesian_gradient(forcefield, positions):
    """Admission check for new element environments; never an MD force method.

    Performed once during preparation. Race forces remain analytic CalcGrad.
    Local checks do not prove accuracy of UFF or cover every thermal geometry.
    """
    flat = np.asarray(positions, dtype=float).ravel().copy()
    gradient = np.asarray(forcefield.CalcGrad(flat.tolist()), dtype=float)
    numerical = np.empty_like(gradient)
    step = 1e-5
    for index in range(len(flat)):
        plus, minus = flat.copy(), flat.copy()
        plus[index] += step; minus[index] -= step
        numerical[index] = (forcefield.CalcEnergy(plus.tolist())-forcefield.CalcEnergy(minus.tolist()))/(2*step)
    differences = np.abs(gradient-numerical)
    tolerance = 1e-3+1e-4*np.maximum(np.abs(gradient), np.abs(numerical))
    if not np.isfinite(numerical).all() or not np.isfinite(gradient).all() or np.any(differences > tolerance):
        worst = int(np.argmax(differences/tolerance))
        raise MoleculePreparationError(
            "当前 UFF 版本对该构象或化学环境的能量与受力不一致，不能开始比赛。"
            f"原子 {worst//3+1} 的 {'xyz'[worst%3]} 分量差 {differences[worst]:.6g} kcal/mol/Å；请调整结构。")
    return {"method": "complete 3N Cartesian central finite difference, once before MD",
            "passed": True, "coordinates_checked": len(flat), "step_angstrom": step,
            "absolute_tolerance_kcal_mol_A": 1e-3, "relative_tolerance": 1e-4,
            "maximum_error_kcal_mol_A": float(differences.max()),
            "qualification": "Local numerical consistency; not chemical accuracy or all-conformation coverage."}


def _read_single_molecule(text: str, filename: str) -> Chem.Mol:
    extension = Path(filename).suffix.lower()
    if extension not in {".mol", ".sdf"}:
        raise MoleculePreparationError("请导入 .mol 或 .sdf 分子文件。")
    if not isinstance(text, str) or not text.strip():
        raise MoleculePreparationError("分子文件为空。")
    if len(text) > 1_000_000:
        raise MoleculePreparationError("文件过大；当前原型每辆车最多支持 256 个原子（含氢）。")

    # Preserve leading blank lines: MOL's title and comment are positional.
    text = text.removeprefix("\ufeff").replace("\r\n", "\n")
    records = [record for record in re.split(r"^\$\$\$\$[ \t]*$", text, flags=re.MULTILINE) if record.strip()]
    if len(records) != 1:
        raise MoleculePreparationError("文件包含多个分子记录。请单独导出一辆赛车的 MOL/SDF 文件。")
    if len(re.findall(r"^M\s+END\s*$", records[0], flags=re.MULTILINE)) != 1:
        raise MoleculePreparationError("分子文件格式不完整，或包含多个结构；请重新导出 MOL/SDF。")

    try:
        with rdBase.BlockLogs():
            molecule = Chem.MolFromMolBlock(records[0], sanitize=False, removeHs=False, strictParsing=True)
        if molecule is None:
            raise MoleculePreparationError("无法读取分子结构，请检查 MOL/SDF 的原子、键和坐标记录。")
        with rdBase.BlockLogs():
            Chem.SanitizeMol(molecule)
    except MoleculePreparationError:
        raise
    except Exception as exc:
        raise MoleculePreparationError(f"化学结构检查失败：请检查价态、键级和芳香环。({exc})") from exc

    if molecule.GetNumAtoms() == 0:
        raise MoleculePreparationError("分子中没有原子。")
    if molecule.GetNumAtoms() > MAX_ATOMS:
        raise MoleculePreparationError("当前原型每辆车最多支持 256 个原子（含氢）。")
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise MoleculePreparationError("赛车必须是一个相互连接的分子；请移除盐、溶剂或未连接的片段。")
    if any(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()):
        raise MoleculePreparationError("当前固定电荷模型尚不支持自由基；请使用闭壳层分子。")
    if any(bond.GetBondType() not in {
        Chem.BondType.SINGLE, Chem.BondType.DOUBLE,
        Chem.BondType.TRIPLE, Chem.BondType.AROMATIC,
    } for bond in molecule.GetBonds()):
        raise MoleculePreparationError("当前原型支持单键、双键、三键和芳香键，尚不支持配位键或查询键。")

    unknown_elements = sorted({atom.GetSymbol() for atom in molecule.GetAtoms()} - EFFECTIVE_POLARIZABILITY.keys())
    if unknown_elements:
        raise MoleculePreparationError("当前驱动模型缺少这些元素的有效极化率参数：" + ", ".join(unknown_elements) + "。")
    if molecule.GetNumConformers() and not np.isfinite(molecule.GetConformer().GetPositions()).all():
        raise MoleculePreparationError("分子坐标包含无效数值，请重新导出结构。")
    return molecule


def prepare_molecule(text: str, filename: str, *, start_pose=None) -> PreparedMolecule:
    """Validate MOL/SDF, add H, prepare 3D geometry and build a single UFF.

    Atom IDs and bond endpoints are 1-based, positions are angstrom, masses
    are atomic mass units, charges are elementary-charge units.  The caller
    must apply force = -CalcGrad(flat_positions); the array order is the
    same as ``atoms`` and must remain stable throughout the simulation.
    The returned force field is mutable and should belong to one race.
    """
    molecule = _read_single_molecule(text, filename)
    title = molecule.GetProp("_Name").strip() if molecule.HasProp("_Name") else ""
    name = (title or Path(filename).stem or "Imported molecule")[:120]
    notes = [
        "分子内作用：UFF；固定电荷：Gasteiger。",
        "有效极化率和针尖驱动使用游戏近似。",
    ]
    # Isotopes remain in the editor/source file, but are outside the race model.
    # Clear labels on this private preparation copy before GetMass, embedding,
    # COM placement and thermal velocities, so all stages use the same masses.
    isotope_count = sum(atom.GetIsotope() != 0 for atom in molecule.GetAtoms())
    if isotope_count:
        for atom in molecule.GetAtoms():
            atom.SetIsotope(0)
        notes.append(f"比赛不区分同位素：已在准备副本中忽略 {isotope_count} 个同位素标签，统一使用元素常规质量。")
    has_3d = bool(molecule.GetNumConformers() and molecule.GetConformer().Is3D())
    original_atom_count = molecule.GetNumAtoms()
    try:
        molecule = Chem.AddHs(molecule, addCoords=has_3d)
    except Exception as exc:
        raise MoleculePreparationError("无法为此结构补齐氢原子；请检查价态和形式电荷。") from exc
    if molecule.GetNumAtoms() > MAX_ATOMS:
        raise MoleculePreparationError(f"补氢后有 {molecule.GetNumAtoms()} 个原子，超过当前 256 个原子的限制。")
    added_hydrogens = molecule.GetNumAtoms() - original_atom_count
    if added_hydrogens:
        notes.append(f"自动补入 {added_hydrogens} 个氢原子。")

    with rdBase.BlockLogs():
        has_uff = AllChem.UFFHasAllMoleculeParams(molecule)
    if not has_uff:
        raise MoleculePreparationError("UFF 缺少此分子中某些原子的参数。请调整结构；游戏不会静默切换力场。")
    if any(atom.GetSymbol() == "P" and atom.GetDegree() == 3 for atom in molecule.GetAtoms()):
        raise MoleculePreparationError(
            "当前 RDKit UFF 版本的三配位磷能量与梯度不一致，暂不支持这类赛车。"
            "四配位磷结构仍可经过赋参与松弛检查后使用。")

    if not has_3d:
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = EMBED_SEED
        parameters.numThreads = 1
        with rdBase.BlockLogs():
            status = AllChem.EmbedMolecule(molecule, parameters)
        if status != 0:
            # The random-coordinate fallback changes the embedding procedure,
            # never the molecular graph or force field, and uses the same seed.
            parameters.useRandomCoords = True
            with rdBase.BlockLogs():
                status = AllChem.EmbedMolecule(molecule, parameters)
        if status != 0:
            raise MoleculePreparationError("无法生成稳定的三维初始构象。请尝试导入已有三维坐标的 MOL 文件。")
        notes.append("从二维结构生成了三维构象（ETKDGv3，固定随机种子）。")
    else:
        notes.append("读取导入的三维构象，做 UFF 几何松弛。")

    try:
        with rdBase.BlockLogs():
            forcefield = AllChem.UFFGetMoleculeForceField(molecule, vdwThresh=100.0, ignoreInterfragInteractions=False)
            forcefield.Initialize()
            result = forcefield.Minimize(maxIts=300, forceTol=0.001, energyTol=1e-7)
        if result != 0:
            raise MoleculePreparationError("初始构象在 300 步 UFF 优化内未收敛；请整理三维结构后再试。")
        positions = np.array(molecule.GetConformer().GetPositions(), dtype=np.float64, copy=True)
        if not np.isfinite(positions).all():
            raise MoleculePreparationError("几何优化产生了无效坐标，请检查初始结构。")
        masses = np.array([atom.GetMass() for atom in molecule.GetAtoms()], dtype=float)
        if not np.isfinite(masses).all() or np.any(masses <= 0):
            raise MoleculePreparationError("某些原子质量无效；请检查元素和同位素信息。")
        positions, pose_report = place_start_pose(positions, masses,
            [atom.GetSymbol() for atom in molecule.GetAtoms()],
            PRESERVE_START_POSE if start_pose is None else start_pose, minimum_z=STARTING_MIN_Z)
        conformer = molecule.GetConformer()
        for index, position in enumerate(positions):
            conformer.SetAtomPosition(index, tuple(float(value) for value in position))
        # Rebuild after translation so the conformer and energy object agree.
        forcefield = AllChem.UFFGetMoleculeForceField(molecule, vdwThresh=100.0, ignoreInterfragInteractions=False)
        forcefield.Initialize()
        flat_coordinates = tuple(float(value) for value in positions.ravel())
        energy = float(forcefield.CalcEnergy(flat_coordinates))
        gradient = np.asarray(forcefield.CalcGrad(flat_coordinates), dtype=float)
        if not np.isfinite(energy) or not np.isfinite(gradient).all():
            raise MoleculePreparationError("此初始结构的 UFF 能量或力无效，不能开始模拟。")
    except MoleculePreparationError:
        raise
    except Exception as exc:
        raise MoleculePreparationError(f"构建 UFF 力场失败：{exc}") from exc

    try:
        with rdBase.BlockLogs():
            AllChem.ComputeGasteigerCharges(molecule, nIter=12, throwOnParamFailure=True)
        charges = np.array([float(atom.GetProp("_GasteigerCharge")) for atom in molecule.GetAtoms()])
    except Exception as exc:
        raise MoleculePreparationError("Gasteiger 电荷模型不能为此结构赋参；请使用支持的普通闭壳层分子。") from exc
    if not np.isfinite(charges).all():
        raise MoleculePreparationError("电荷计算产生了无效数值，不能开始模拟。")
    formal_charge = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    if not np.isclose(charges.sum(), formal_charge, atol=1e-5):
        raise MoleculePreparationError("分子部分电荷的总和与形式电荷不一致，不能开始模拟。")
    if formal_charge:
        notes.append(f"此分子净电荷为 {formal_charge:+d} e，电场响应与中性分子不同。")

    has_hetero = any(atom.GetSymbol() not in ("C", "H") for atom in molecule.GetAtoms())
    if has_hetero and (abs(energy) > 1e8 or np.linalg.norm(gradient.reshape(-1, 3), axis=1).max() > 1e5):
        raise MoleculePreparationError(
            "UFF 优化后仍有异常大的能量或受力，当前构象或化学环境暂不能用于比赛；请调整结构。")
    forcefield_validation = (_check_uff_cartesian_gradient(forcefield, positions) if has_hetero
                            else {"method": "C/H analytic UFF path retained; separate release numerical checks",
                                  "qualification": "No new per-import finite-difference gate for the established C/H path."})

    vdw_parameters = []
    for atom in molecule.GetAtoms():
        with rdBase.BlockLogs():
            parameters = AllChem.GetUFFVdWParams(molecule, atom.GetIdx(), atom.GetIdx())
        if parameters is None or not np.isfinite(parameters).all() or min(parameters) <= 0:
            raise MoleculePreparationError(f"UFF 无法为原子 {atom.GetIdx()+1}（{atom.GetSymbol()}）提供有效的范德华参数。")
        vdw_parameters.append(tuple(float(value) for value in parameters))

    atoms = [
        {
            "id": index + 1,
            "element": atom.GetSymbol(),
            "mass": float(masses[index]),
            "charge": float(charges[index]),
            "polarizability": EFFECTIVE_POLARIZABILITY[atom.GetSymbol()],
            "uff_vdw_rmin_angstrom": vdw_parameters[index][0],
            "uff_vdw_well_depth_kcal_mol": vdw_parameters[index][1],
        }
        for index, atom in enumerate(molecule.GetAtoms())
    ]
    bonds = [
        {"a": bond.GetBeginAtomIdx() + 1, "b": bond.GetEndAtomIdx() + 1, "order": float(bond.GetBondTypeAsDouble())}
        for bond in molecule.GetBonds()
    ]
    return PreparedMolecule(name, atoms, bonds, positions, forcefield, notes, molecule, forcefield_validation, pose_report)
