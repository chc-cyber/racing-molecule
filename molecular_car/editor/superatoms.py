"""SUP display groups carried by explicit file row / stable editor ID maps.

These groups never replace chemistry with pseudo-atoms. All members and bonds
remain in the molecule; collapsing changes only its two-dimensional depiction.
"""
from uuid import uuid4
from rdkit import Chem


def read_groups(molecule, index_to_atom_id):
    result = []
    for group in Chem.GetMolSubstanceGroups(molecule):
        if group.GetProp("TYPE") != "SUP":
            raise ValueError("当前支持原生 Superatom 折叠基团；其他 S-group 尚不能无损编辑。")
        members = [index_to_atom_id[index] for index in group.GetAtoms()]
        if not members:
            raise ValueError("Superatom 基团没有真实内部原子。")
        points = []
        for point in group.GetAttachPoints():
            points.append({"atomId": index_to_atom_id[point.aIdx],
                           "attachmentId": str(point.id),
                           **({"leavingAtomId": index_to_atom_id[point.lvIdx]} if point.lvIdx >= 0 else {})})
        result.append({"id": "g-" + str(uuid4()), "type": "SUP",
                       "name": group.GetProp("LABEL") if group.HasProp("LABEL") else "Group",
                       "atomIds": members, "attachmentPoints": points,
                       "expanded": group.HasProp("ESTATE") and group.GetProp("ESTATE") == "E"})
    return result


def write_groups(molecule, project):
    indices = {atom["id"]: index for index, atom in enumerate(project["atoms"])}
    for record in project.get("metadata", {}).get("ketcherGroups", []):
        members = record.get("atomIds", [])
        if record.get("type") != "SUP" or not members or len(set(members)) != len(members) or any(i not in indices for i in members):
            raise ValueError("折叠基团的真实原子或编号映射无效。")
        group = Chem.CreateMolSubstanceGroup(molecule, "SUP")
        group.SetProp("LABEL", record.get("name", "Group"))
        if record.get("expanded", False):
            group.SetProp("ESTATE", "E")
        member_indices = {indices[i] for i in members}
        for i in members:
            group.AddAtomWithIdx(indices[i])
        for bond in molecule.GetBonds():
            if (bond.GetBeginAtomIdx() in member_indices) != (bond.GetEndAtomIdx() in member_indices):
                group.AddBondWithIdx(bond.GetIdx())
        for point in record.get("attachmentPoints", []):
            if point.get("atomId") not in members:
                raise ValueError("折叠基团的连接点不属于其真实原子。")
            leaving = point.get("leavingAtomId")
            if leaving is not None and leaving not in indices:
                raise ValueError("折叠基团的离去原子引用无效。")
            group.AddAttachPoint(indices[point["atomId"]], indices[leaving] if leaving is not None else -1,
                                 str(point.get("attachmentId", "1")))
    return molecule
