from __future__ import annotations

import math

from .constants import (
    LIMIT_ENERGY,
    XS_TYPE_DUMMY,
    XS_TYPE_H,
    XS_TYPE_SIZE,
    xs_hbond,
    xs_is_heavy,
    xs_is_hydrophobic,
    xs_radius,
)
from .model import Atom, Molecule

THRESHOLD = 8
PRECI = 1000
SZ = THRESHOLD * PRECI
TERM_WEIGHTS = (-0.035579, -0.005156, 0.840245, -0.035069, -0.587439)


def sqr(value):
    return value * value


class EnergyCalculator:
    def __init__(self, rad_scale=1.0):
        self.rad_scale = float(rad_scale)
        self.precalculated = {}

    def _surface_distance(self, t1, t2, distance):
        return distance - (xs_radius(t1) + xs_radius(t2)) * self.rad_scale

    def gauss1(self, t1, t2, d):
        if t1 in {XS_TYPE_H, XS_TYPE_DUMMY} or t2 in {XS_TYPE_H, XS_TYPE_DUMMY}:
            return 0.0
        return math.exp(-sqr(d * 2.0))

    def gauss2(self, t1, t2, d):
        if t1 in {XS_TYPE_H, XS_TYPE_DUMMY} or t2 in {XS_TYPE_H, XS_TYPE_DUMMY}:
            return 0.0
        return math.exp(-sqr((d - 3.0) * 0.5))

    def repulsion(self, t1, t2, d):
        return 0.0 if d > 0 else sqr(d)

    def hydrophobic(self, t1, t2, d):
        if t1 in {XS_TYPE_H, XS_TYPE_DUMMY} or t2 in {XS_TYPE_H, XS_TYPE_DUMMY}:
            return 0.0
        if not xs_is_hydrophobic(t1) or not xs_is_hydrophobic(t2):
            return 0.0
        return 0.0 if d >= 1.5 else (1.0 if d <= 0.5 else 1.5 - d)

    def hydrogen_bond(self, t1, t2, d):
        if t1 in {XS_TYPE_H, XS_TYPE_DUMMY} or t2 in {XS_TYPE_H, XS_TYPE_DUMMY}:
            return 0.0
        if not xs_hbond(t1, t2):
            return 0.0
        return 0.0 if d >= 0.0 else (1.0 if d <= -0.7 else d * -1.428571)

    def _weighted(self, t1, t2, d):
        return (
            TERM_WEIGHTS[0] * self.gauss1(t1, t2, d)
            + TERM_WEIGHTS[1] * self.gauss2(t1, t2, d)
            + TERM_WEIGHTS[2] * self.repulsion(t1, t2, d)
            + TERM_WEIGHTS[3] * self.hydrophobic(t1, t2, d)
            + TERM_WEIGHTS[4] * self.hydrogen_bond(t1, t2, d)
        )

    def _precalculated_atom_energy(self, atom1, atom2):
        t1, t2 = atom1.xs_type, atom2.xs_type
        distance = atom1.distance_to(atom2)
        if distance > THRESHOLD:
            return 0.0
        idx = int(distance * PRECI)
        if idx >= SZ:
            return 0.0
        key = (t1, t2, idx)
        if key not in self.precalculated:
            d = self._surface_distance(t1, t2, idx / float(PRECI))
            self.precalculated[key] = self._weighted(t1, t2, d)
        return self.precalculated[key]

    def get_energy_strict(self, atom1, atom2):
        distance = atom1.distance_to(atom2)
        if distance > THRESHOLD:
            return 0.0
        d = self._surface_distance(atom1.xs_type, atom2.xs_type, distance)
        return self._weighted(atom1.xs_type, atom2.xs_type, d)

    def get_energy(self, atom_or_mol, atom_or_mol2):
        if isinstance(atom_or_mol, Atom) and isinstance(atom_or_mol2, Atom):
            return self._precalculated_atom_energy(atom_or_mol, atom_or_mol2)
        if isinstance(atom_or_mol, Atom) and isinstance(atom_or_mol2, Molecule):
            total = 0.0
            for atom in atom_or_mol2.atoms:
                if atom.xs_type != XS_TYPE_H:
                    total += self.get_energy(atom_or_mol, atom)
                    if total >= LIMIT_ENERGY:
                        return LIMIT_ENERGY
            return total
        if isinstance(atom_or_mol, Molecule) and isinstance(atom_or_mol2, Molecule):
            total = 0.0
            for atom in atom_or_mol.atoms:
                if atom.xs_type != XS_TYPE_H:
                    total += self.get_energy(atom, atom_or_mol2)
                    if total >= LIMIT_ENERGY:
                        return LIMIT_ENERGY
            return total
        raise TypeError("unsupported energy operands")

    def calc_intra_energy(self, ligand):
        total = 0.0
        distances = ligand.graph_distances()
        for i, atom_i in enumerate(ligand.atoms):
            if not xs_is_heavy(atom_i.xs_type):
                continue
            for atom_j in ligand.atoms[i + 1 :]:
                if not xs_is_heavy(atom_j.xs_type):
                    continue
                if distances[atom_i.id][atom_j.id] != -1 and distances[atom_i.id][atom_j.id] < 4:
                    continue
                total += self.get_energy(atom_i, atom_j)
                if total >= LIMIT_ENERGY:
                    return LIMIT_ENERGY
        return total

    def affinity(self, ligand, receptor):
        score = self.get_energy(ligand, receptor)
        return score / (1.0 + 0.05846 * ligand.nrots())
