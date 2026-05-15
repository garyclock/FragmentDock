from pathlib import Path
import tempfile
import unittest

from restretto.config import DockingConfiguration, SearchGrid
from restretto.energy import EnergyCalculator
import restretto.docking_full_rotation as docking_full_rotation_module
from restretto.fragment_grid_6d import build_fragment_6d_grid
from restretto.geometry import Point3d, Vector3d
from restretto.model import Atom, Molecule
from restretto.rotation import (
    choose_better_rotation,
    make_initial_rotations,
    nearest_rotation_bin,
    rotation_distance,
)
from restretto.constants import XS_TYPE_C_H


DATA = Path("references/restretto/testdata")


class FullRotationTests(unittest.TestCase):
    def test_initial_rotations_are_stable_and_nonempty(self):
        rotations = make_initial_rotations()

        self.assertEqual(len(rotations), 60)
        self.assertEqual(rotations[0].as_tuple(), (0.0, 0.0, 0.0))
        self.assertEqual([rot.as_tuple() for rot in rotations], [rot.as_tuple() for rot in make_initial_rotations()])

    def test_nearest_rotation_bin_uses_energy_tie_break(self):
        bins = [Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 2.0)]

        idx = nearest_rotation_bin(Vector3d(0.0, 0.0, 1.0), bins, energies=[5.0, -2.0])

        self.assertEqual(idx, 1)
        self.assertAlmostEqual(rotation_distance(Vector3d(0.0, 0.0, 1.0), bins[0]), 1.0)
        self.assertEqual(choose_better_rotation((0, 1.0, 5.0), (1, 1.0, -2.0)), (1, 1.0, -2.0))

    def test_fragment_6d_grid_scores_translation_and_rotation(self):
        fragment = Molecule([Atom(0, Vector3d(0.0, 0.0, 0.0), XS_TYPE_C_H)], "frag", "C")
        receptor = Molecule([Atom(0, Vector3d(3.8, 0.0, 0.0), XS_TYPE_C_H)], "rec", "C")
        rotations = [Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 1.0)]

        grid = build_fragment_6d_grid(
            fragment,
            receptor,
            translations=[Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0)],
            rotations=rotations,
            energy_calculator=EnergyCalculator(),
        )

        self.assertEqual(grid.score_pose(0, 0), grid.query_nearest(Point3d(0.0, 0.0, 0.0), rotations[0]).energy)
        self.assertLessEqual(grid.query_nearest(Point3d(0.1, 0.0, 0.0), Vector3d(0.0, 0.0, 0.5)).rotation_idx, 1)

    def test_full_rotation_considers_candidates_beyond_the_initial_prefix(self):
        tmpdir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(tmpdir, ignore_errors=True))
        rotangs = tmpdir / "rotangs.csv"
        rotangs.write_text(
            "\n".join(
                [
                    "0.0, 0.0, 0.0",
                    "0.0, 0.0, 3.141592653589793",
                    "0.0, 0.0, 1.5707963267948966",
                ]
            ),
            encoding="utf-8",
        )

        receptor = Molecule(
            [
                Atom(0, Vector3d(0.0, 1.0, 0.0), XS_TYPE_C_H),
            ],
            "rec",
            "C",
        )
        ligand = Molecule(
            [
                Atom(0, Vector3d(-1.0, 0.0, 0.0), XS_TYPE_C_H),
                Atom(1, Vector3d(1.0, 0.0, 0.0), XS_TYPE_C_H),
            ],
            "lig",
            "C",
        )

        conf = DockingConfiguration(
            grid=SearchGrid(
                center=Point3d(0.0, 0.0, 0.0),
                outer_width=Point3d(0.0, 0.0, 0.0),
                inner_width=Point3d(0.0, 0.0, 0.0),
                search_pitch=Point3d(1.0, 1.0, 1.0),
                score_pitch=Point3d(1.0, 1.0, 1.0),
            ),
            receptor_file=str(DATA / "2HU4_A_r.pdbqt"),
            output_file=str(tmpdir / "out.sdf"),
            rotangs_file=str(rotangs),
            poses_per_lig=1,
            poses_per_lig_before_opt=2,
            no_local_opt=True,
        )

        class FakeEnergyCalculator(object):
            def __init__(self, rad_scale=1.0):
                self.rad_scale = rad_scale

            def affinity(self, ligand_mol, receptor_mol):
                for atom in ligand_mol.atoms:
                    if atom.y > 0.9:
                        return -1.0
                return 0.0

        original_ec = docking_full_rotation_module.EnergyCalculator
        docking_full_rotation_module.EnergyCalculator = FakeEnergyCalculator
        try:
            scored = docking_full_rotation_module.dock_full_rotation(conf, receptor, [ligand])
        finally:
            docking_full_rotation_module.EnergyCalculator = original_ec

        self.assertEqual(len(scored), 1)
        best_pose, _ = scored[0]
        closest = min(best_pose.atoms, key=lambda atom: atom.distance_to(receptor.atoms[0]))
        self.assertLess(abs(closest.y - 1.0), 0.2)


if __name__ == "__main__":
    unittest.main()
