import math
import unittest

from restretto.constants import XS_TYPE_C_H, XS_TYPE_H, XS_TYPE_N_D, XS_TYPE_O_A
from restretto.geometry import Vector3d
from restretto.model import Atom, Bond, Molecule


class GeometryModelTests(unittest.TestCase):
    def assertTupleAlmostEqual(self, left, right, places=6):
        self.assertEqual(len(left), len(right))
        for lval, rval in zip(left, right):
            self.assertAlmostEqual(lval, rval, places=places)

    def test_vector_axis_rotation_around_z(self):
        vec = Vector3d(1.0, 0.0, 0.0)

        rotated = vec.axis_rotated(Vector3d(0.0, 0.0, 1.0), math.pi / 2)

        self.assertAlmostEqual(rotated.x, 0.0, places=6)
        self.assertAlmostEqual(rotated.y, 1.0, places=6)
        self.assertAlmostEqual(rotated.z, 0.0, places=6)


    def test_molecule_center_radius_graph_distance_and_hydrogen_deletion(self):
        mol = Molecule(
            atoms=[
                Atom(0, Vector3d(0.0, 0.0, 0.0), XS_TYPE_C_H),
                Atom(1, Vector3d(2.0, 0.0, 0.0), XS_TYPE_N_D),
                Atom(2, Vector3d(4.0, 0.0, 0.0), XS_TYPE_O_A),
                Atom(3, Vector3d(9.0, 0.0, 0.0), XS_TYPE_H),
            ],
            title="chain",
            smiles="CNO",
        )
        mol.append_bond(Bond(0, 1, True))
        mol.append_bond(Bond(1, 2, True))
        mol.append_bond(Bond(2, 3, False))

        self.assertTupleAlmostEqual(mol.center().as_tuple(), (2.0, 0.0, 0.0))
        self.assertAlmostEqual(mol.radius(), 2.0)
        self.assertEqual(mol.graph_distances()[0][2], 2)

        mol.delete_hydrogens()

        self.assertEqual(mol.size, 3)
        self.assertEqual([atom.id for atom in mol.atoms], [0, 1, 2])
        self.assertEqual(mol.graph_distances()[0][2], 2)


if __name__ == "__main__":
    unittest.main()
