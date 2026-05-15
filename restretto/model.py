from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass
from typing import List

from .constants import XS_TYPE_DUMMY, XS_TYPE_H, xs_is_heavy
from .geometry import Vector3d


@dataclass
class Atom(Vector3d):
    id: int = 0
    xs_type: int = 0

    def __init__(self, atom_id, pos, xs_type):
        super().__init__(pos.x, pos.y, pos.z)
        self.id = int(atom_id)
        self.xs_type = int(xs_type)

    def translated(self, vec):
        return Atom(self.id, Vector3d(self.x + vec.x, self.y + vec.y, self.z + vec.z), self.xs_type)

    def rotated_atom(self, theta, phi, psi):
        return Atom(self.id, self.rotated(theta, phi, psi), self.xs_type)

    def axis_rotated_atom(self, axis, theta):
        return Atom(self.id, self.axis_rotated(axis, theta), self.xs_type)


@dataclass
class Bond:
    atom_id1: int
    atom_id2: int
    is_rotor: bool = False
    order: int = 1


class Molecule:
    def __init__(self, atoms=None, title="", smiles=""):
        self.atoms: List[Atom] = list(atoms or [])
        self.bonds: List[Bond] = []
        self.title = title
        self.smiles = smiles
        self.identifier = "%s,%s" % (title, smiles)
        self.intra_energy = 0.0
        self.bond_ids = []
        self._rebuild_bond_ids()

    @property
    def size(self):
        return len(self.atoms)

    def _rebuild_bond_ids(self):
        max_id = max((atom.id for atom in self.atoms), default=-1)
        self.bond_ids = [[] for _ in range(max_id + 1)]
        for idx, bond in enumerate(self.bonds):
            needed = max(bond.atom_id1, bond.atom_id2) + 1
            while len(self.bond_ids) < needed:
                self.bond_ids.append([])
            self.bond_ids[bond.atom_id1].append(idx)
            self.bond_ids[bond.atom_id2].append(idx)

    def append_atom(self, atom):
        self.atoms.append(atom)
        while len(self.bond_ids) <= atom.id:
            self.bond_ids.append([])

    def append_bond(self, bond):
        idx = len(self.bonds)
        self.bonds.append(bond)
        while len(self.bond_ids) <= max(bond.atom_id1, bond.atom_id2):
            self.bond_ids.append([])
        self.bond_ids[bond.atom_id1].append(idx)
        self.bond_ids[bond.atom_id2].append(idx)

    def copy(self):
        clone = Molecule([Atom(atom.id, atom, atom.xs_type) for atom in self.atoms], self.title, self.smiles)
        for bond in self.bonds:
            clone.append_bond(Bond(bond.atom_id1, bond.atom_id2, bond.is_rotor, getattr(bond, "order", 1)))
        clone.intra_energy = self.intra_energy
        return clone

    def translate(self, vec):
        self.atoms = [atom.translated(vec) for atom in self.atoms]

    def rotate(self, theta, phi, psi):
        self.atoms = [atom.rotated_atom(theta, phi, psi) for atom in self.atoms]

    def rotated_copy(self, rotation):
        clone = self.copy()
        center = clone.center()
        clone.translate(Vector3d(-center.x, -center.y, -center.z))
        clone.rotate(rotation.x, rotation.y, rotation.z)
        clone.translate(center)
        return clone

    def axis_rotate(self, axis, theta):
        self.atoms = [atom.axis_rotated_atom(axis, theta) for atom in self.atoms]

    def center(self):
        heavy = [atom for atom in self.atoms if atom.xs_type not in {XS_TYPE_H, XS_TYPE_DUMMY}]
        if not heavy:
            return Vector3d()
        total = Vector3d()
        for atom in heavy:
            total = total + atom
        return total / len(heavy)

    def radius(self):
        center = self.center()
        return max(((atom - center).abs() for atom in self.atoms if atom.xs_type not in {XS_TYPE_H, XS_TYPE_DUMMY}), default=0.0)

    def heavy_num(self):
        return sum(1 for atom in self.atoms if atom.xs_type != XS_TYPE_H)

    def graph_distances(self):
        n = max((atom.id for atom in self.atoms), default=-1) + 1
        dist = [[-1 for _ in range(n)] for _ in range(n)]
        for atom in self.atoms:
            start = atom.id
            dist[start][start] = 0
            queue = deque([start])
            while queue:
                point = queue.popleft()
                if point >= len(self.bond_ids):
                    continue
                for bond_idx in self.bond_ids[point]:
                    bond = self.bonds[bond_idx]
                    other = bond.atom_id1 + bond.atom_id2 - point
                    if dist[start][other] == -1:
                        dist[start][other] = dist[start][point] + 1
                        queue.append(other)
        return dist

    def delete_hydrogens(self):
        old_to_new = {}
        new_atoms = []
        for atom in self.atoms:
            if atom.xs_type == XS_TYPE_H:
                continue
            new_id = len(new_atoms)
            old_to_new[atom.id] = new_id
            new_atoms.append(Atom(new_id, atom, atom.xs_type))
        new_bonds = []
        for bond in self.bonds:
            if bond.atom_id1 in old_to_new and bond.atom_id2 in old_to_new:
                new_bonds.append(Bond(old_to_new[bond.atom_id1], old_to_new[bond.atom_id2], bond.is_rotor, getattr(bond, "order", 1)))
        self.atoms = new_atoms
        self.bonds = new_bonds
        self._rebuild_bond_ids()

    def renumbering(self, new_size, labels):
        new_atoms = [None for _ in range(int(new_size))]
        old_id_to_index = {atom.id: idx for idx, atom in enumerate(self.atoms)}
        for idx, label in enumerate(labels):
            if label == 0:
                continue
            new_idx = int(label) - 1
            atom = self.atoms[idx]
            new_atoms[new_idx] = Atom(new_idx, atom, atom.xs_type)
        new_bonds = []
        for bond in self.bonds:
            if bond.atom_id1 not in old_id_to_index or bond.atom_id2 not in old_id_to_index:
                continue
            left_label = labels[old_id_to_index[bond.atom_id1]]
            right_label = labels[old_id_to_index[bond.atom_id2]]
            if left_label == 0 or right_label == 0:
                continue
            new_bonds.append(Bond(int(left_label) - 1, int(right_label) - 1, bond.is_rotor, getattr(bond, "order", 1)))
        self.atoms = new_atoms
        self.bonds = new_bonds
        self._rebuild_bond_ids()

    def nrots(self):
        n = max((atom.id for atom in self.atoms), default=-1) + 1
        adj_heavy = [0 for _ in range(n)]
        atom_by_id = {atom.id: atom for atom in self.atoms}
        for atom in self.atoms:
            for bond_idx in self.bond_ids[atom.id]:
                bond = self.bonds[bond_idx]
                other_id = bond.atom_id1 + bond.atom_id2 - atom.id
                other = atom_by_id[other_id]
                if xs_is_heavy(other.xs_type):
                    adj_heavy[atom.id] += 1
        ret = 0.0
        for atom in self.atoms:
            if atom.xs_type == XS_TYPE_H:
                continue
            count = 0
            for bond_idx in self.bond_ids[atom.id]:
                bond = self.bonds[bond_idx]
                other_id = bond.atom_id1 + bond.atom_id2 - atom.id
                other = atom_by_id[other_id]
                if bond.is_rotor and xs_is_heavy(other.xs_type) and adj_heavy[other.id] > 1:
                    count += 1
            ret += count * 0.5
        return ret

    def calc_rmsd(self, other):
        if self.size != other.size:
            return math.inf
        sd = 0.0
        count = 0
        for a, b in zip(self.atoms, other.atoms):
            if a.xs_type in {XS_TYPE_H, XS_TYPE_DUMMY} or b.xs_type in {XS_TYPE_H, XS_TYPE_DUMMY}:
                continue
            sd += (a - b).norm()
            count += 1
        return math.sqrt(sd / count) if count else 0.0
