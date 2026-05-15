from __future__ import annotations

from .constants import XS_TYPE_DUMMY, XS_TYPE_H
from .model import Atom, Bond, Molecule


class _UnionFind:
    def __init__(self, size):
        self.parent = list(range(size))

    def find(self, value):
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra

    def same(self, a, b):
        return self.find(a) == self.find(b)

    def size(self, value):
        root = self.find(value)
        return sum(1 for i in range(len(self.parent)) if self.find(i) == root)

    def sets(self):
        groups = {}
        for i in range(len(self.parent)):
            groups.setdefault(self.find(i), []).append(i)
        return list(groups.values())


def _ring_detector(atom_num, bonds):
    edges = [[] for _ in range(atom_num)]
    for bond in bonds:
        edges[bond.atom_id1].append(bond.atom_id2)
        edges[bond.atom_id2].append(bond.atom_id1)

    rings = []
    for start in range(atom_num):
        done = [False] * atom_num
        route = []
        ring = []

        def dfs(now):
            nonlocal ring
            if done[now]:
                if now == start and len(route) > 2:
                    ring = list(route)
                return
            done[now] = True
            route.append(now)
            for nxt in edges[now]:
                dfs(nxt)
            route.pop()

        dfs(start)
        if ring:
            rings.append(ring)
    return rings


def _extract_substructure(mol, atom_ids):
    id_map = {old: new for new, old in enumerate(atom_ids)}
    sub = Molecule([Atom(id_map[old], mol.atoms[old], mol.atoms[old].xs_type) for old in atom_ids], mol.title, mol.smiles)
    for bond in mol.bonds:
        if bond.atom_id1 in id_map and bond.atom_id2 in id_map:
            sub.append_bond(Bond(id_map[bond.atom_id1], id_map[bond.atom_id2], bond.is_rotor, getattr(bond, "order", 1)))
    return sub


def _bond_rotate(mol, bond_id, theta):
    if bond_id >= len(mol.bonds):
        raise IndexError("invalid bond id")
    uf = _UnionFind(mol.size)
    for idx, bond in enumerate(mol.bonds):
        if idx != bond_id:
            uf.union(bond.atom_id1, bond.atom_id2)
    sets = uf.sets()
    if sets and len(sets[0]) == mol.size:
        return mol
    bond = mol.bonds[bond_id]
    axis = mol.atoms[bond.atom_id2] - mol.atoms[bond.atom_id1]
    new_mol = mol.copy()
    anchor = mol.atoms[bond.atom_id1]
    new_mol.translate(-anchor)
    for atom_id in sets[0]:
        new_mol.atoms[atom_id] = new_mol.atoms[atom_id].axis_rotated_atom(axis, theta)
    new_mol.translate(mol.center() - new_mol.center())
    return new_mol


def _is_mergeable(mol, atomids_a, atomids_b):
    united_ids = list(atomids_a) + list(atomids_b)
    united = _extract_substructure(mol, united_ids)
    prev_a = _extract_substructure(mol, atomids_a)
    prev_b = _extract_substructure(mol, atomids_b)
    prev_rings = len(_ring_detector(prev_a.size, prev_a.bonds)) + len(_ring_detector(prev_b.size, prev_b.bonds))
    united_rings = len(_ring_detector(united.size, united.bonds))
    if prev_rings != united_rings:
        return False
    for idx, bond in enumerate(united.bonds):
        if not bond.is_rotor:
            continue
        rotated = _bond_rotate(united, idx, 1.0)
        if united.calc_rmsd(rotated) >= 1e-5:
            return False
    return True


def _find_set_containing(sets, atom_id):
    for group in sets:
        if atom_id in group:
            return group
    return []


def decompose_molecule(molecule, max_ring_size=-1, merge_solitary=True):
    atoms = molecule.atoms
    bonds = molecule.bonds
    uf = _UnionFind(molecule.size)

    for bond in bonds:
        if not bond.is_rotor:
            uf.union(bond.atom_id1, bond.atom_id2)

    for ring in _ring_detector(len(atoms), bonds):
        if max_ring_size != -1 and len(ring) > max_ring_size:
            continue
        for atom_id in ring:
            uf.union(ring[0], atom_id)

    heavy_adj_count = [0 for _ in atoms]
    for bond in bonds:
        a, b = bond.atom_id1, bond.atom_id2
        if atoms[a].xs_type != XS_TYPE_H:
            heavy_adj_count[b] += 1
        if atoms[b].xs_type != XS_TYPE_H:
            heavy_adj_count[a] += 1

    done = [False for _ in atoms]
    for bond in bonds:
        a, b = bond.atom_id1, bond.atom_id2
        if atoms[a].xs_type == XS_TYPE_H or atoms[b].xs_type == XS_TYPE_H:
            continue
        for _ in range(2):
            if (
                (uf.size(a) == 1 and heavy_adj_count[a] <= 2 and uf.size(b) > 1 and not done[b])
                or (uf.size(a) == 1 and uf.size(b) == 1 and heavy_adj_count[a] <= 2 and heavy_adj_count[b] <= 2)
            ):
                uf.union(a, b)
                done[a] = True
            a, b = b, a

    for bond in bonds:
        if not merge_solitary:
            break
        a, b = bond.atom_id1, bond.atom_id2
        if atoms[a].xs_type == XS_TYPE_H or atoms[b].xs_type == XS_TYPE_H:
            continue
        if uf.same(a, b):
            continue
        sets = uf.sets()
        atomids_a = _find_set_containing(sets, a)
        atomids_b = _find_set_containing(sets, b)
        if _is_mergeable(molecule, atomids_a, atomids_b):
            uf.union(a, b)

    for bond in bonds:
        a, b = bond.atom_id1, bond.atom_id2
        if atoms[a].xs_type == XS_TYPE_H or atoms[b].xs_type == XS_TYPE_H:
            uf.union(a, b)

    id_sets = uf.sets()
    set_id = [0 for _ in atoms]
    for idx, group in enumerate(id_sets):
        for atom_id in group:
            set_id[atom_id] = idx

    dummies = [[] for _ in id_sets]
    bonds_in_fragments = [[] for _ in id_sets]
    for bond in bonds:
        a, b = bond.atom_id1, bond.atom_id2
        if atoms[a].xs_type == XS_TYPE_H or atoms[b].xs_type == XS_TYPE_H:
            continue
        if set_id[a] != set_id[b]:
            dummies[set_id[a]].append(Atom(b, atoms[b], XS_TYPE_DUMMY))
            dummies[set_id[b]].append(Atom(a, atoms[a], XS_TYPE_DUMMY))
            bonds_in_fragments[set_id[a]].append(Bond(a, b, bond.is_rotor, getattr(bond, "order", 1)))
            bonds_in_fragments[set_id[b]].append(Bond(a, b, bond.is_rotor, getattr(bond, "order", 1)))
        else:
            bonds_in_fragments[set_id[a]].append(Bond(a, b, bond.is_rotor, getattr(bond, "order", 1)))

    fragments = []
    for idx, group in enumerate(id_sets):
        fragment_atoms = []
        sorted_dummies = sorted(dummies[idx], key=lambda atom: atom.id)
        dummy_i = 0
        for atom_id in group:
            if atoms[atom_id].xs_type == XS_TYPE_H:
                continue
            while dummy_i < len(sorted_dummies) and sorted_dummies[dummy_i].id < atoms[atom_id].id:
                fragment_atoms.append(sorted_dummies[dummy_i])
                dummy_i += 1
            fragment_atoms.append(Atom(atoms[atom_id].id, atoms[atom_id], atoms[atom_id].xs_type))
        while dummy_i < len(sorted_dummies):
            fragment_atoms.append(sorted_dummies[dummy_i])
            dummy_i += 1

        frag = Molecule(fragment_atoms, "%s_frag_%d" % (molecule.title, idx), molecule.smiles)
        for bond in bonds_in_fragments[idx]:
            frag.append_bond(bond)
        fragments.append(frag)
    return fragments
