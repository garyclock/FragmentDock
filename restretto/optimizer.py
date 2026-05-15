from __future__ import annotations

import math

from .constants import LIMIT_ENERGY, XS_TYPE_H
from .geometry import Vector3d


class _CRand:
    RAND_MAX = 2147483647

    def __init__(self, seed=0):
        seed = 1 if int(seed) == 0 else int(seed)
        self._state = [0 for _ in range(344)]
        self._state[0] = seed
        for i in range(1, 31):
            self._state[i] = (16807 * self._state[i - 1]) % self.RAND_MAX
        for i in range(31, 34):
            self._state[i] = self._state[i - 31]
        for i in range(34, 344):
            self._state[i] = (self._state[i - 31] + self._state[i - 3]) & 0xFFFFFFFF
        self._idx = 344

    def rand(self):
        value = (self._state[self._idx - 31] + self._state[self._idx - 3]) & 0xFFFFFFFF
        self._state.append(value)
        self._idx += 1
        return value >> 1

    def randf(self, lower, upper):
        return (self.rand() / self.RAND_MAX) * (upper - lower) + lower


class OptimizerGrid:
    def __init__(self, atom_grids, max_rmsd=1e10):
        self.atom_grids = list(atom_grids)
        self.max_rmsd = float(max_rmsd)

    def calc_inter_energy(self, mol):
        total = 0.0
        for atom in mol.atoms:
            if atom.xs_type == XS_TYPE_H:
                continue
            total += self.atom_grids[atom.xs_type].get_inter_energy(atom)
            if total >= LIMIT_ENERGY:
                return LIMIT_ENERGY
        return total

    def calc_total_energy(self, mol):
        return self.calc_inter_energy(mol) + mol.intra_energy

    def optimize(self, mol):
        rng = _CRand(0)
        opt = self.calc_total_energy(mol)
        nearest_num = 200
        trans_step = 0.5
        rotate_step = math.pi / 30.0
        initial_mol = mol.copy()

        while True:
            next_val = 1e10
            next_mol = mol.copy()
            center = mol.center()
            for _ in range(nearest_num):
                dv = Vector3d(
                    rng.randf(-trans_step, trans_step),
                    rng.randf(-trans_step, trans_step),
                    rng.randf(-trans_step, trans_step),
                )
                theta = rng.randf(-rotate_step, rotate_step)
                phi = rng.randf(-rotate_step, rotate_step)
                psi = rng.randf(-rotate_step, rotate_step)

                tmp = mol.copy()
                tmp.translate(-center)
                tmp.rotate(theta, phi, psi)
                tmp.translate(center + dv)

                if initial_mol.calc_rmsd(tmp) > self.max_rmsd:
                    continue

                val = self.calc_total_energy(tmp)
                if val < next_val:
                    next_val = val
                    next_mol = tmp

            if next_val < opt:
                opt = next_val
                mol.atoms = next_mol.atoms
                mol.bonds = next_mol.bonds
                mol.bond_ids = next_mol.bond_ids
            else:
                break

        return opt
