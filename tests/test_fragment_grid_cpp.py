import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch

from restretto.constants import LIMIT_ENERGY, XS_TYPE_C_H
from restretto.fragment_grid import FragmentInterEnergyGrid, FragmentInterEnergyGridContainer, make_distance_grid
from restretto.fragment_grid_6d import FragmentInterEnergyGrid6D, build_fragment_rotation_cache, nearest_fragment_rotation_bin
from restretto.geometry import Point3d, Vector3d
from restretto.grid import InterEnergyGrid
from restretto.model import Atom, Bond, Molecule
from restretto.config import DockingConfiguration, ReuseStrategy


class FragmentGridCppTests(unittest.TestCase):
    def test_single_atom_fragment_grid_matches_atom_grid_values(self):
        atom_grid = InterEnergyGrid(Point3d(0, 0, 0), Point3d(1, 1, 1), Point3d(3, 3, 3), LIMIT_ENERGY)
        atom_grid.set_inter_energy(1, 1, 1, -2.0)
        atom_grids = [atom_grid for _ in range(21)]
        receptor = Molecule([Atom(0, Vector3d(3, 0, 0), XS_TYPE_C_H)], "rec", "C")
        distance_grid = make_distance_grid(atom_grid.center, atom_grid.pitch, atom_grid.num, receptor)
        fragment = Molecule([Atom(0, Vector3d(0, 0, 0), XS_TYPE_C_H)], "frag", "C")
        fragment.idx = 7

        fg = FragmentInterEnergyGrid(fragment, [Vector3d(0, 0, 0)], atom_grids, distance_grid)

        self.assertEqual(fg.frag_idx, 7)
        self.assertAlmostEqual(fg.grid.get_inter_energy(1, 1, 1), -2.0)
        self.assertEqual(fg.grid.get_inter_energy(0, 0, 0), LIMIT_ENERGY)

    def test_single_atom_6d_fragment_grid_matches_atom_grid_values(self):
        atom_grid = InterEnergyGrid(Point3d(0, 0, 0), Point3d(1, 1, 1), Point3d(3, 3, 3), LIMIT_ENERGY)
        atom_grid.set_inter_energy(1, 1, 1, -2.0)
        atom_grids = [atom_grid for _ in range(21)]
        receptor = Molecule([Atom(0, Vector3d(3, 0, 0), XS_TYPE_C_H)], "rec", "C")
        distance_grid = make_distance_grid(atom_grid.center, atom_grid.pitch, atom_grid.num, receptor)
        fragment = Molecule([Atom(0, Vector3d(0, 0, 0), XS_TYPE_C_H)], "frag", "C")
        fragment.idx = 7

        fg = FragmentInterEnergyGrid6D(fragment, [Vector3d(0, 0, 0), Vector3d(0, 0, 1.0)], atom_grids, distance_grid)

        self.assertEqual(fg.frag_idx, 7)
        self.assertEqual(fg.rotation_count, 1)
        self.assertAlmostEqual(fg.get_inter_energy(0, 1, 1, 1), -2.0)
        self.assertEqual(fg.get_inter_energy(0, 0, 0, 0), LIMIT_ENERGY)

    def test_6d_fragment_grid_preserves_rotation_specific_slices(self):
        atom_grid = InterEnergyGrid(Point3d(0, 0, 0), Point3d(1, 1, 1), Point3d(5, 5, 5), LIMIT_ENERGY)
        atom_grid.set_inter_energy(2, 2, 2, -1.0)
        atom_grid.set_inter_energy(3, 2, 2, -10.0)
        atom_grid.set_inter_energy(2, 3, 2, -100.0)
        atom_grids = [atom_grid for _ in range(21)]
        receptor = Molecule([Atom(0, Vector3d(3, 0, 0), XS_TYPE_C_H)], "rec", "C")
        distance_grid = make_distance_grid(atom_grid.center, atom_grid.pitch, atom_grid.num, receptor)
        fragment = Molecule(
            [
                Atom(0, Vector3d(0, 0, 0), XS_TYPE_C_H),
                Atom(1, Vector3d(1, 0, 0), XS_TYPE_C_H),
            ],
            "frag",
            "CC",
        )
        fragment.idx = 3

        fg = FragmentInterEnergyGrid6D(
            fragment,
            [Vector3d(0, 0, 0), Vector3d(0, 0, 1.5707963267948966)],
            atom_grids,
            distance_grid,
        )

        self.assertAlmostEqual(fg.get_inter_energy(0, 2, 2, 2), -11.0)
        self.assertAlmostEqual(fg.get_inter_energy(1, 2, 2, 2), -101.0)

    def test_nearest_fragment_rotation_bin_uses_geometry_and_tie_breaks_low_index(self):
        fragment = Molecule(
            [
                Atom(0, Vector3d(0, 0, 0), XS_TYPE_C_H),
                Atom(1, Vector3d(1, 0, 0), XS_TYPE_C_H),
            ],
            "frag",
            "CC",
        )
        occurrence = fragment.copy()
        occurrence.rotate(0, 0, 1.5707963267948966)

        idx = nearest_fragment_rotation_bin(
            fragment,
            occurrence,
            [Vector3d(0, 0, 0), Vector3d(0, 0, 1.5707963267948966)],
        )
        tie_idx = nearest_fragment_rotation_bin(fragment, fragment.copy(), [Vector3d(0, 0, 0), Vector3d(0, 0, 0)])

        self.assertEqual(idx, 1)
        self.assertEqual(tie_idx, 0)

    def test_cached_nearest_fragment_rotation_bin_matches_geometry_lookup(self):
        fragment = Molecule(
            [
                Atom(0, Vector3d(0, 0, 0), XS_TYPE_C_H),
                Atom(1, Vector3d(1, 0, 0), XS_TYPE_C_H),
                Atom(2, Vector3d(0, 1, 0), XS_TYPE_C_H),
            ],
            "frag",
            "CCC",
        )
        occurrence = fragment.copy()
        occurrence.rotate(0.2, 0.1, 0.9)
        rotations = [
            Vector3d(0, 0, 0),
            Vector3d(0.2, 0.1, 0.9),
            Vector3d(0.0, 0.0, 1.5707963267948966),
        ]

        cache = build_fragment_rotation_cache(fragment, rotations)

        self.assertEqual(nearest_fragment_rotation_bin(cache, occurrence), 1)

    def test_fragment_grid_container_online_reuses_registered_grid(self):
        grid_a = FragmentInterEnergyGrid.__new__(FragmentInterEnergyGrid)
        grid_a.frag_idx = 1
        grid_b = FragmentInterEnergyGrid.__new__(FragmentInterEnergyGrid)
        grid_b.frag_idx = 2
        container = FragmentInterEnergyGridContainer(1)

        container.insert(grid_a)
        self.assertTrue(container.is_registered(1))
        self.assertIs(container.get(1), grid_a)
        container.next()
        container.insert(grid_b)

        self.assertFalse(container.is_registered(1))
        self.assertTrue(container.is_registered(2))

    def test_fragment_grid_container_online_get_refreshes_lru_slot(self):
        grid_a = FragmentInterEnergyGrid.__new__(FragmentInterEnergyGrid)
        grid_a.frag_idx = 1
        grid_b = FragmentInterEnergyGrid.__new__(FragmentInterEnergyGrid)
        grid_b.frag_idx = 2
        grid_c = FragmentInterEnergyGrid.__new__(FragmentInterEnergyGrid)
        grid_c.frag_idx = 3
        container = FragmentInterEnergyGridContainer(2)

        container.insert(grid_a)
        container.next()
        container.insert(grid_b)
        container.next()
        self.assertIs(container.get(1), grid_a)
        container.next()
        container.insert(grid_c)

        self.assertTrue(container.is_registered(1))
        self.assertFalse(container.is_registered(2))
        self.assertTrue(container.is_registered(3))

    def test_fragment_grid_container_offline_uses_scheduled_slot_only(self):
        grid_a = FragmentInterEnergyGrid.__new__(FragmentInterEnergyGrid)
        grid_a.frag_idx = 1
        grid_b = FragmentInterEnergyGrid.__new__(FragmentInterEnergyGrid)
        grid_b.frag_idx = 2
        container = FragmentInterEnergyGridContainer(2, indices_to_save=[1, 0])

        container.insert(grid_a)
        self.assertTrue(container.is_registered(1))
        container.next()

        self.assertFalse(container.is_registered(1))
        container.insert(grid_b)
        self.assertIs(container.get(2), grid_b)

    def test_offline_schedule_matches_cpp_left_back_ssp_for_repeated_fragment(self):
        from restretto.fragment_grid import build_offline_schedule

        schedule = build_offline_schedule([[(0, 2), (1, 3), (0, 2)]], cache_size=1, fragment_count=2)

        self.assertEqual(schedule, [0, 0, 0])

    def test_fragment_grid_docking_scores_single_atom_like_cpp_fraggrid_path(self):
        from restretto.cli import _dock_with_fragment_grids

        tmpdir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(tmpdir, ignore_errors=True))
        grid_dir = tmpdir / "grid"
        grid_dir.mkdir()
        atom_grid = InterEnergyGrid(Point3d(0, 0, 0), Point3d(1, 1, 1), Point3d(1, 1, 1), 0.0)
        atom_grid.set_inter_energy(0, 0, 0, -5.0)
        from restretto.constants import XS_STRINGS

        for name in XS_STRINGS:
            atom_grid.write_grid(grid_dir / ("%s.grid" % name))

        conf = DockingConfiguration()
        conf.grid.center = Point3d(0, 0, 0)
        conf.grid.inner_width = Point3d(0, 0, 0)
        conf.grid.outer_width = Point3d(0, 0, 0)
        conf.grid.search_pitch = Point3d(1, 1, 1)
        conf.grid.score_pitch = Point3d(1, 1, 1)
        conf.grid_folder = str(grid_dir)
        conf.reuse_grid = ReuseStrategy.NONE
        conf.output_score_threshold = 0.0
        conf.no_local_opt = True
        receptor = Molecule([Atom(0, Vector3d(10, 10, 10), XS_TYPE_C_H)], "rec", "C")
        ligand = Molecule([Atom(0, Vector3d(0, 0, 0), XS_TYPE_C_H)], "lig", "C")

        scored = _dock_with_fragment_grids(tmpdir / "testgrid.in", conf, receptor, [ligand])

        self.assertEqual(len(scored), 1)
        self.assertAlmostEqual(scored[0][1], -5.0)

    def test_full_rotation_docking_keeps_rotation_specific_fragment_grid_scores(self):
        from restretto.cli import _dock_with_fragment_grids, _dock_with_fragment_grids_6d
        from restretto.constants import XS_STRINGS

        tmpdir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(tmpdir, ignore_errors=True))
        grid_dir = tmpdir / "grid"
        grid_dir.mkdir()
        atom_grid = InterEnergyGrid(Point3d(0, 0, 0), Point3d(1, 1, 1), Point3d(7, 5, 5), LIMIT_ENERGY)
        atom_grid.set_inter_energy(2, 2, 2, -1.0)
        atom_grid.set_inter_energy(3, 2, 2, -20.0)
        atom_grid.set_inter_energy(4, 2, 2, -10.0)
        atom_grid.set_inter_energy(5, 2, 2, -30.0)
        atom_grid.set_inter_energy(3, 1, 2, -100.0)
        atom_grid.set_inter_energy(3, 3, 2, -1000.0)
        for name in XS_STRINGS:
            atom_grid.write_grid(grid_dir / ("%s.grid" % name))

        conf = DockingConfiguration()
        conf.grid.center = Point3d(0, 0, 0)
        conf.grid.inner_width = Point3d(2, 0, 0)
        conf.grid.outer_width = Point3d(6, 4, 4)
        conf.grid.search_pitch = Point3d(1, 1, 1)
        conf.grid.score_pitch = Point3d(1, 1, 1)
        conf.grid_folder = str(grid_dir)
        conf.reuse_grid = ReuseStrategy.NONE
        conf.output_score_threshold = -40.0
        conf.no_local_opt = True
        conf.poses_per_lig_before_opt = 1
        receptor = Molecule([Atom(0, Vector3d(3, 0, 0), XS_TYPE_C_H)], "rec", "C")
        ligand = Molecule(
            [
                Atom(0, Vector3d(0, 0, 0), XS_TYPE_C_H),
                Atom(1, Vector3d(1, 0, 0), XS_TYPE_C_H),
            ],
            "lig",
            "CC",
        )
        ligand.append_bond(Bond(0, 1, False))
        rotations = [Vector3d(0, 0, 0), Vector3d(0, 0, 1.5707963267948966)]

        with patch("restretto.cli.make_initial_rotations", return_value=[Vector3d(0, 0, 0)]), patch(
            "restretto.cli.make_rotations_60", return_value=rotations
        ), patch("restretto.cli.nearest_fragment_rotation_bin", return_value=0):
            default_scored = _dock_with_fragment_grids(tmpdir / "testgrid.in", conf, receptor, [ligand])
            full_scored = _dock_with_fragment_grids_6d(tmpdir / "testgrid.in", conf, receptor, [ligand])

        self.assertLess(default_scored[0][0].center().x, full_scored[0][0].center().x)

    def test_fragment_grid_docking_places_all_output_poses_in_search_box_coordinates(self):
        from restretto.cli import main
        from restretto.io import read_molecules

        tmpdir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(tmpdir, ignore_errors=True))
        data = Path("references/restretto/testdata")
        grid_dir = tmpdir / "grid"
        output = tmpdir / "docked.sdf"
        config = tmpdir / "testgrid.in"
        config.write_text(
            "\n".join(
                [
                    "INNERBOX 2, 2, 2",
                    "OUTERBOX 2, 2, 2",
                    "BOX_CENTER 0.3826, 81.705, 109.195",
                    "SEARCH_PITCH 1, 1, 1",
                    "SCORING_PITCH 1, 1, 1",
                    "MEMORY_SIZE 10",
                    f"RECEPTOR {data / '2HU4_A_r.pdbqt'}",
                    f"LIGAND {data / 'G39.mol2'}",
                    f"OUTPUT {output}",
                    f"GRID_FOLDER {grid_dir}",
                    "NO_LOCAL_OPT true",
                    "POSES_PER_LIG 2",
                    "POSES_PER_LIG_BEFORE_OPT 2",
                    "OUTPUT_SCORE_THRESHOLD 100",
                ]
            ),
            encoding="utf-8",
        )

        self.assertEqual(main(["atomgrid-gen", str(config)]), 0)
        self.assertEqual(main(["conformer-docking", str(config)]), 0)

        mols = read_molecules(output)
        self.assertGreaterEqual(len(mols), 1)
        for mol in mols:
            center = mol.center()
            self.assertLess(abs(center.x - 0.3826), 3.0)
            self.assertLess(abs(center.y - 81.705), 3.0)
            self.assertLess(abs(center.z - 109.195), 3.0)

    def test_select_output_poses_deduplicates_rmsd_similar_poses(self):
        from restretto.cli import _select_output_poses

        mol_a = Molecule(
            [
                Atom(0, Vector3d(0.0, 0.0, 0.0), XS_TYPE_C_H),
                Atom(1, Vector3d(1.5, 0.0, 0.0), XS_TYPE_C_H),
            ],
            "lig",
            "CC",
        )
        mol_b = mol_a.copy()
        mol_c = mol_a.copy()
        mol_c.translate(Vector3d(5.0, 0.0, 0.0))

        selected = _select_output_poses(
            [(-3.0, mol_a), (-2.5, mol_b), (-2.0, mol_c)],
            pose_min_rmsd=0.5,
            max_poses=2,
        )

        self.assertEqual(len(selected), 2)
        self.assertIs(selected[0][1], mol_a)
        self.assertIs(selected[1][1], mol_c)

    def test_normalized_fragment_applies_cpp_reference_triplet_rotation(self):
        from restretto.cli import _normalized_fragment

        fragment = Molecule(
            [
                Atom(0, Vector3d(1, 1, 1), XS_TYPE_C_H),
                Atom(1, Vector3d(1, 2, 1), XS_TYPE_C_H),
                Atom(2, Vector3d(1, 1, 2), XS_TYPE_C_H),
            ],
            "bent",
            "CCC",
        )
        fragment.append_bond(__import__("restretto.model", fromlist=["Bond"]).Bond(0, 1, False))
        fragment.append_bond(__import__("restretto.model", fromlist=["Bond"]).Bond(0, 2, False))

        normalized = _normalized_fragment(fragment)
        tri1 = normalized.atoms[1] - normalized.atoms[0]
        tri2 = normalized.atoms[2] - normalized.atoms[0]

        self.assertAlmostEqual(normalized.center().abs(), 0.0, places=6)
        self.assertAlmostEqual(tri1.y, 0.0, places=6)
        self.assertAlmostEqual(tri1.z, 0.0, places=6)
        self.assertAlmostEqual(tri2.z, 0.0, places=6)


if __name__ == "__main__":
    unittest.main()
