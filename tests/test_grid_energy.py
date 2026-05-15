from pathlib import Path
import tempfile
import unittest

from restretto.constants import LIMIT_ENERGY, XS_TYPE_C_H, XS_TYPE_N_D, XS_TYPE_O_A
from restretto.energy import EnergyCalculator
from restretto.geometry import Point3d, Vector3d
from restretto.grid import InterEnergyGrid
from restretto.model import Atom, Bond, Molecule


class GridEnergyTests(unittest.TestCase):
    def assertTupleAlmostEqual(self, left, right, places=6):
        self.assertEqual(len(left), len(right))
        for lval, rval in zip(left, right):
            self.assertAlmostEqual(lval, rval, places=places)

    def test_grid_roundtrip_and_bounds(self):
        grid = InterEnergyGrid(Point3d(0.0, 0.0, 0.0), Point3d(0.5, 0.5, 0.5), Point3d(3, 3, 3), 0.0)
        grid.set_inter_energy(1, 1, 1, -2.5)
        grid.set_inter_energy(2, 1, 0, 4.25)

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "small.grid"
            grid.write_grid(path)
            loaded = InterEnergyGrid.parse_grid(path)

        self.assertAlmostEqual(loaded.get_inter_energy(1, 1, 1), -2.5)
        self.assertAlmostEqual(loaded.get_inter_energy(2, 1, 0), 4.25)
        self.assertEqual(loaded.get_inter_energy(-1, 0, 0), LIMIT_ENERGY)
        self.assertTupleAlmostEqual(loaded.convert(2, 1, 0).as_tuple(), (0.5, 0.0, -0.5))

    def test_grid_exposes_numpy_values_for_cpp_grid_slices(self):
        grid = InterEnergyGrid(Point3d(0, 0, 0), Point3d(1, 1, 1), Point3d(2, 2, 2), 0.0)
        grid.set_inter_energy(1, 0, 1, -3.0)

        values = grid.values3d()

        self.assertEqual(values.shape, (2, 2, 2))
        self.assertAlmostEqual(float(values[1, 0, 1]), -3.0)


    def test_energy_strict_and_precalculated_are_close(self):
        carbon = Atom(0, Vector3d(0.0, 0.0, 0.0), XS_TYPE_C_H)
        donor = Atom(1, Vector3d(2.8, 0.0, 0.0), XS_TYPE_N_D)
        acceptor = Atom(2, Vector3d(6.1, 0.0, 0.0), XS_TYPE_O_A)
        ec = EnergyCalculator()

        self.assertAlmostEqual(ec.get_energy(carbon, donor), ec.get_energy_strict(carbon, donor), delta=0.002)
        self.assertGreater(ec.hydrogen_bond(XS_TYPE_N_D, XS_TYPE_O_A, -0.5), 0)
        self.assertLess(ec.get_energy(donor, acceptor), 0)


    def test_intra_energy_skips_atoms_closer_than_four_bonds(self):
        mol = Molecule(
            [
                Atom(0, Vector3d(0.0, 0.0, 0.0), XS_TYPE_C_H),
                Atom(1, Vector3d(1.0, 0.0, 0.0), XS_TYPE_C_H),
                Atom(2, Vector3d(2.0, 0.0, 0.0), XS_TYPE_C_H),
                Atom(3, Vector3d(3.0, 0.0, 0.0), XS_TYPE_C_H),
                Atom(4, Vector3d(2.5, 0.0, 0.0), XS_TYPE_O_A),
            ],
            "chain",
            "CCCCC",
        )
        for idx in range(4):
            mol.append_bond(Bond(idx, idx + 1, True))

        ec = EnergyCalculator()

        self.assertAlmostEqual(ec.calc_intra_energy(mol), ec.get_energy(mol.atoms[0], mol.atoms[4]))


if __name__ == "__main__":
    unittest.main()
