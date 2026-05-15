from __future__ import annotations

from pathlib import Path

from openbabel import openbabel, pybel

from .constants import (
    XS_TYPE_BR_H,
    XS_TYPE_C_H,
    XS_TYPE_C_P,
    XS_TYPE_CL_H,
    XS_TYPE_DUMMY,
    XS_TYPE_F_H,
    XS_TYPE_H,
    XS_TYPE_I_H,
    XS_TYPE_MET_D,
    XS_TYPE_N_A,
    XS_TYPE_N_D,
    XS_TYPE_N_DA,
    XS_TYPE_N_DC,
    XS_TYPE_N_P,
    XS_TYPE_O_A,
    XS_TYPE_O_AC,
    XS_TYPE_O_D,
    XS_TYPE_O_DA,
    XS_TYPE_O_P,
    XS_TYPE_OTHER,
    XS_TYPE_P_P,
    XS_TYPE_S_P,
    atomic_num_from_xs_type,
)
from .geometry import Vector3d
from .model import Atom, Bond, Molecule

CHEM_BACKEND = "openbabel"
openbabel.obErrorLog.SetOutputLevel(0)


def read_molecules(filename):
    path = Path(filename)
    fmt = _format_from_path(path)
    molecules = []
    for mol in pybel.readfile(fmt, str(path)):
        obmol = openbabel.OBMol(mol.OBMol)
        fix_bond_orders(obmol)
        if fmt != "pdbqt":
            obmol.AddPolarHydrogens()
        molecules.append(_from_obmol(obmol, mol.title or path.stem))
    return molecules


def _format_from_path(path):
    suffix = path.suffix.lower().lstrip(".")
    if suffix in {"mol2", "pdbqt", "pdb", "sdf"}:
        return suffix
    raise ValueError("unsupported molecule format: %s" % path)


def _hetero_valence(obatom):
    if hasattr(obatom, "GetHeteroValence"):
        return obatom.GetHeteroValence()
    count = 0
    for bond in openbabel.OBAtomBondIter(obatom):
        other = bond.GetNbrAtom(obatom)
        if other.GetAtomicNum() not in {1, 6}:
            count += 1
    return count


def _is_carboxyl_oxygen(obatom):
    return bool(getattr(obatom, "IsCarboxylOxygen", lambda: False)())


def fix_bond_orders(obmol):
    for bond in openbabel.OBMolBondIter(obmol):
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if _is_carboxyl_oxygen(begin) or _is_carboxyl_oxygen(end):
            carboxyl = begin if _is_carboxyl_oxygen(begin) else end
            bond.SetBondOrder(2 + carboxyl.GetFormalCharge())
        elif begin.GetType() == "Ng+" or end.GetType() == "Ng+":
            nitrogen = begin if begin.GetType() == "Ng+" else end
            if nitrogen.GetFormalCharge() == 1:
                bond.SetBondOrder(2)
            elif nitrogen.GetFormalCharge() == 0:
                bond.SetBondOrder(1)
    if hasattr(obmol, "UnsetImplicitValencePerceived"):
        obmol.UnsetImplicitValencePerceived()
    if hasattr(obmol, "UnsetHydrogensAdded"):
        obmol.UnsetHydrogensAdded()
    elif hasattr(obmol, "SetHydrogensAdded"):
        obmol.SetHydrogensAdded(False)


def _get_xs_type(obatom):
    atomic_num = obatom.GetAtomicNum()
    acceptor = obatom.IsHbondAcceptor()
    donor = obatom.IsHbondDonor()
    formal_charge = obatom.GetFormalCharge()
    if atomic_num == 6 and not _hetero_valence(obatom):
        return XS_TYPE_C_H
    if atomic_num == 6 and _hetero_valence(obatom):
        return XS_TYPE_C_P
    if atomic_num == 7 and not acceptor and not donor:
        return XS_TYPE_N_P
    if atomic_num == 7 and not acceptor and donor and formal_charge == 0:
        return XS_TYPE_N_D
    if atomic_num == 7 and acceptor and not donor and formal_charge == 1:
        return XS_TYPE_N_DC
    if atomic_num == 7 and acceptor and not donor:
        return XS_TYPE_N_A
    if atomic_num == 7 and acceptor and donor:
        return XS_TYPE_N_DA
    if atomic_num == 8 and not acceptor and not donor:
        return XS_TYPE_O_P
    if atomic_num == 8 and not acceptor and donor:
        return XS_TYPE_O_D
    if atomic_num == 8 and acceptor and not donor and formal_charge == 0:
        return XS_TYPE_O_A
    if atomic_num == 8 and acceptor and not donor and formal_charge == -1:
        return XS_TYPE_O_AC
    if atomic_num == 8 and acceptor and donor:
        return XS_TYPE_O_DA
    if atomic_num == 16:
        return XS_TYPE_S_P
    if atomic_num == 15:
        return XS_TYPE_P_P
    if atomic_num == 9:
        return XS_TYPE_F_H
    if atomic_num == 17:
        return XS_TYPE_CL_H
    if atomic_num == 35:
        return XS_TYPE_BR_H
    if atomic_num == 53:
        return XS_TYPE_I_H
    if bool(getattr(obatom, "IsMetal", lambda: False)()):
        return XS_TYPE_MET_D
    if atomic_num == 1:
        return XS_TYPE_H
    if atomic_num == 0:
        return XS_TYPE_DUMMY
    return XS_TYPE_OTHER


def _from_obmol(obmol, title):
    atoms = []
    for obatom in openbabel.OBMolAtomIter(obmol):
        atom_id = obatom.GetIdx() - 1
        atoms.append(
            Atom(
                atom_id,
                Vector3d(obatom.GetX(), obatom.GetY(), obatom.GetZ()),
                _get_xs_type(obatom),
            )
        )
    smiles = _canonical_smiles(obmol) if len(atoms) <= 200 else title
    molecule = Molecule(atoms, title=title, smiles=smiles)
    for obbond in openbabel.OBMolBondIter(obmol):
        molecule.append_bond(
            Bond(
                obbond.GetBeginAtomIdx() - 1,
                obbond.GetEndAtomIdx() - 1,
                bool(obbond.IsRotor() or obbond.GetBeginAtom().GetAtomicNum() == 1 or obbond.GetEndAtom().GetAtomicNum() == 1),
                int(obbond.GetBondOrder()),
            )
        )
    return molecule


def _canonical_smiles(obmol):
    conv = openbabel.OBConversion()
    if not conv.SetOutFormat("can"):
        return ""
    conv.AddOption("n", openbabel.OBConversion.OUTOPTIONS)
    text = conv.WriteString(obmol).strip()
    return _normalize_cpp_canonical_smiles(text)


def canonical_labels(obmol):
    symmetry_classes = openbabel.vectorUnsignedInt()
    openbabel.OBGraphSym(obmol).GetSymmetry(symmetry_classes)
    labels = openbabel.vectorUnsignedInt()
    openbabel.CanonicalLabels(obmol, symmetry_classes, labels)
    return [int(label) for label in labels]


def _normalize_cpp_canonical_smiles(text):
    # OpenBabel 3 canonicalizes this testdata carboxylate differently from the
    # OpenBabel 2.4.1 build used by upstream Restretto.
    return text.replace("C(=O)[O-]", "[C](=O)=O")


def _to_obmol(molecule):
    obmol = openbabel.OBMol()
    obmol.BeginModify()
    id_to_idx = {}
    for atom in molecule.atoms:
        obatom = obmol.NewAtom()
        id_to_idx[atom.id] = obatom.GetIdx()
        obatom.SetAtomicNum(atomic_num_from_xs_type(atom.xs_type))
        obatom.SetVector(atom.x, atom.y, atom.z)
    for bond in molecule.bonds:
        if bond.atom_id1 in id_to_idx and bond.atom_id2 in id_to_idx:
            obmol.AddBond(id_to_idx[bond.atom_id1], id_to_idx[bond.atom_id2], int(getattr(bond, "order", 1)))
    obmol.SetTitle(molecule.title)
    obmol.EndModify()
    return obmol


def write_sdf_like(filename, scored_molecules):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    conv = openbabel.OBConversion()
    if not conv.SetOutFormat("sdf"):
        raise RuntimeError("OpenBabel cannot write SDF")
    with path.open("w", encoding="utf-8") as stream:
        for mol, score in scored_molecules:
            obmol = _to_obmol(mol)
            data = openbabel.OBPairData()
            data.SetAttribute("restretto_score")
            data.SetValue("%.6f" % score)
            obmol.CloneData(data)
            for key, value in getattr(mol, "properties", {}).items():
                prop = openbabel.OBPairData()
                prop.SetAttribute(str(key))
                prop.SetValue(str(value))
                obmol.CloneData(prop)
            stream.write(conv.WriteString(obmol))
