import unittest

from restretto.constants import XS_TYPE_C_H
from restretto.geometry import Point3d, Vector3d
from restretto.grid import InterEnergyGrid
from restretto.model import Atom, Molecule
from restretto.optimizer import OptimizerGrid
from restretto.optimizer import OptimizerFragmentGrid6D
from restretto.optimizer import _CRand


class OptimizerGridTests(unittest.TestCase):
    def test_c_rand_matches_linux_reference_sequence(self):
        rng = _CRand(0)

        self.assertEqual([rng.rand() for _ in range(3)], [1804289383, 846930886, 1681692777])

    def test_grid_optimizer_reduces_total_energy_and_moves_pose(self):
        grid = InterEnergyGrid(Point3d(0.0, 0.0, 0.0), Point3d(1.0, 1.0, 1.0), Point3d(5, 5, 5), 10.0)
        grid.set_inter_energy(2, 2, 2, -5.0)
        atom_grids = [grid for _ in range(21)]
        mol = Molecule([Atom(0, Vector3d(0.6, 0.0, 0.0), XS_TYPE_C_H)], "lig", "C")
        mol.intra_energy = 0.0

        optimizer = OptimizerGrid(atom_grids)
        before = optimizer.calc_total_energy(mol)
        after = optimizer.optimize(mol)

        self.assertLessEqual(after, before)
        self.assertLessEqual(after, -5.0)
        self.assertLess(abs(mol.atoms[0].x), 0.75)

    def test_fragment_grid_6d_optimizer_uses_fragment_grid_scoring(self):
        class FakeFragmentGrid:
            rotation_count = 1

            def get_inter_energy(self, rotation_idx, x, y, z):
                return -5.0 if (x, y, z) == (2, 2, 2) else 10.0

        score_grid = InterEnergyGrid(Point3d(0.0, 0.0, 0.0), Point3d(1.0, 1.0, 1.0), Point3d(5, 5, 5), 0.0)
        base = Molecule([Atom(0, Vector3d(0.0, 0.0, 0.0), XS_TYPE_C_H)], "lig", "C")
        base.intra_energy = 0.0
        frag_item = {"fragment": base.copy(), "frag_idx": 0}

        optimizer = OptimizerFragmentGrid6D(
            score_grid,
            [FakeFragmentGrid()],
            [base.copy()],
            [frag_item],
            [Vector3d(0.0, 0.0, 0.0)],
        )

        pose = base.copy()
        before = optimizer.calc_total_energy(pose)
        after = optimizer.optimize(pose)

        self.assertEqual(before, -5.0)
        self.assertLessEqual(after, before)


if __name__ == "__main__":
    unittest.main()
