from pathlib import Path
import math
import shutil
import tempfile
import unittest

from restretto.cli import build_parser, main
from restretto.config import DockingConfiguration
from restretto.constants import XS_STRINGS, XS_TYPE_SIZE
from restretto.geometry import Point3d, round_half_up, ceili
from restretto.grid import InterEnergyGrid
from restretto.io import read_molecules
from restretto.rotation import make_rotations_60
from restretto.constants import XS_TYPE_C_P, XS_TYPE_O_AC
from restretto.constants import XS_TYPE_DUMMY
from restretto.decompose import decompose_molecule


DATA = Path("references/restretto/testdata")


class CppParityTests(unittest.TestCase):
    def test_docking_configuration_defaults_match_cpp(self):
        conf = DockingConfiguration()

        self.assertTrue(conf.reorder)
        self.assertEqual(conf.poses_per_lig, 1)
        self.assertEqual(conf.poses_per_lig_before_opt, 2000)
        self.assertAlmostEqual(conf.output_score_threshold, -3.0)
        self.assertAlmostEqual(conf.pose_min_rmsd, 0.5)
        self.assertAlmostEqual(conf.local_max_rmsd, 1e10)
        self.assertAlmostEqual(conf.rad_scale, 0.95)

    def test_cpp_round_and_ceil_semantics(self):
        self.assertEqual(round_half_up(0.49), 0)
        self.assertEqual(round_half_up(0.5), 1)
        self.assertEqual(round_half_up(-0.49), 0)
        self.assertEqual(round_half_up(-0.5), -1)
        self.assertEqual(ceili(3.0 - 1e-7), 3)
        self.assertEqual(ceili(3.0 + 1e-7), 4)

    def test_make_rotations_60_matches_reference_algorithm_shape(self):
        rotations = make_rotations_60()

        self.assertEqual(len(rotations), 60)
        self.assertFalse(any(math.isnan(v) for rot in rotations for v in rot.as_tuple()))
        self.assertEqual(rotations[0].as_tuple(), (-0.0, -0.0, 0.0))

    def test_atomgrid_gen_writes_all_cpp_named_xscore_grids(self):
        tmpdir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(tmpdir, ignore_errors=True))
        config = tmpdir / "testgrid.in"
        grid_dir = tmpdir / "grid"
        config.write_text(
            "\n".join(
                [
                    "INNERBOX 4, 4, 4",
                    "OUTERBOX 4, 4, 4",
                    "BOX_CENTER 0.3826, 81.705, 109.195",
                    "SEARCH_PITCH 1, 1, 1",
                    "SCORING_PITCH 1, 1, 1",
                    "MEMORY_SIZE 10",
                    f"RECEPTOR {DATA / '2HU4_A_r.pdbqt'}",
                    f"LIGAND {DATA / 'G39.mol2'}",
                    f"OUTPUT {tmpdir / 'docked.sdf'}",
                    f"GRID_FOLDER {grid_dir}",
                ]
            ),
            encoding="utf-8",
        )

        self.assertEqual(main(["atomgrid-gen", str(config)]), 0)

        grid_files = sorted(path.name for path in grid_dir.glob("*.grid"))
        self.assertEqual(grid_files, sorted(f"{name}.grid" for name in XS_STRINGS))
        self.assertEqual(len(grid_files), XS_TYPE_SIZE)
        parsed = InterEnergyGrid.parse_grid(grid_dir / "C_H.grid")
        self.assertEqual(parsed.num, Point3d(5, 5, 5))

    def test_vectorized_atomgrid_matches_scalar_energy_calculator(self):
        from restretto.cli import _build_atom_grid
        from restretto.energy import EnergyCalculator
        from restretto.model import Atom, Molecule

        receptor = Molecule(
            [
                Atom(0, Point3d(0.0, 0.0, 0.0), XS_TYPE_C_P),
                Atom(1, Point3d(1.2, 0.0, 0.0), XS_TYPE_O_AC),
            ],
            "rec",
            "",
        )
        center = Point3d(0.0, 0.0, 0.0)
        pitch = Point3d(0.5, 0.5, 0.5)
        num = Point3d(3, 3, 3)
        grid = _build_atom_grid(center, pitch, num, XS_TYPE_C_P, receptor, EnergyCalculator(0.95))
        scalar = EnergyCalculator(0.95)

        for x in range(int(num.x)):
            for y in range(int(num.y)):
                for z in range(int(num.z)):
                    pos = grid.convert(x, y, z)
                    probe = Atom(0, pos, XS_TYPE_C_P)
                    self.assertAlmostEqual(grid.get_inter_energy(x, y, z), scalar.get_energy(probe, receptor), places=5)

    def test_openbabel_reading_adds_polar_hydrogens_and_cpp_atom_types(self):
        ligand = read_molecules(DATA / "G39.mol2")[0]

        self.assertEqual(ligand.size, 23)
        self.assertGreaterEqual(sum(1 for atom in ligand.atoms if atom.xs_type == XS_TYPE_C_P), 5)
        self.assertGreaterEqual(sum(1 for atom in ligand.atoms if atom.xs_type == XS_TYPE_O_AC), 1)

    def test_openbabel_canonical_identifier_matches_cpp_carboxylate_form(self):
        ligand = read_molecules(DATA / "G39.mol2")[0]

        self.assertEqual(ligand.smiles, "CCC(O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)[C](=O)=O)CC")

    def test_decompose_adds_dummy_atoms_on_fragment_boundaries(self):
        ligand = read_molecules(DATA / "G39.mol2")[0]
        ligand.delete_hydrogens()

        fragments = decompose_molecule(ligand)

        self.assertGreater(len(fragments), 1)
        self.assertTrue(any(atom.xs_type == XS_TYPE_DUMMY for fragment in fragments for atom in fragment.atoms))

    def test_sdf_output_uses_cpp_score_property_name(self):
        from restretto.io import write_sdf_like
        from restretto.model import Atom, Bond, Molecule
        from restretto.constants import XS_TYPE_C_H

        tmpdir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(tmpdir, ignore_errors=True))
        mol = Molecule(
            [
                Atom(0, Point3d(0.0, 0.0, 0.0), XS_TYPE_C_H),
                Atom(1, Point3d(1.2, 0.0, 0.0), XS_TYPE_C_H),
            ],
            "ethene",
            "C=C",
        )
        mol.append_bond(Bond(0, 1, False, 2))
        output = tmpdir / "out.sdf"

        write_sdf_like(output, [(mol, -1.25)])

        text = output.read_text(encoding="utf-8")
        self.assertIn(">  <restretto_score>", text)
        self.assertNotIn("REstrettoPyScore", text)
        self.assertIn("  1  2  2", text)

    def test_score_only_output_uses_cpp_labels(self):
        import contextlib
        import io

        tmpdir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(tmpdir, ignore_errors=True))
        config = tmpdir / "testgrid.in"
        config.write_text(
            "\n".join(
                [
                    "INNERBOX 4, 4, 4",
                    "OUTERBOX 4, 4, 4",
                    "BOX_CENTER 0, 0, 0",
                    "SEARCH_PITCH 1, 1, 1",
                    "SCORING_PITCH 1, 1, 1",
                    "MEMORY_SIZE 10",
                    f"RECEPTOR {DATA / '2HU4_A_r.pdbqt'}",
                    f"LIGAND {DATA / 'G39.mol2'}",
                    f"OUTPUT {tmpdir / 'out.sdf'}",
                    f"GRID_FOLDER {tmpdir / 'grid'}",
                ]
            ),
            encoding="utf-8",
        )

        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            self.assertEqual(main(["score-only", str(config)]), 0)

        text = stdout.getvalue()
        self.assertIn("Title: ", text)
        self.assertIn("Affinity: ", text)

    def test_conformer_docking_requires_atom_grids_in_default_mode(self):
        tmpdir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(tmpdir, ignore_errors=True))
        config = tmpdir / "testgrid.in"
        config.write_text(
            "\n".join(
                [
                    "INNERBOX 2, 2, 2",
                    "OUTERBOX 2, 2, 2",
                    "BOX_CENTER 0, 0, 0",
                    "SEARCH_PITCH 1, 1, 1",
                    "SCORING_PITCH 1, 1, 1",
                    "MEMORY_SIZE 10",
                    f"RECEPTOR {DATA / '2HU4_A_r.pdbqt'}",
                    f"LIGAND {DATA / 'G39.mol2'}",
                    f"OUTPUT {tmpdir / 'out.sdf'}",
                    f"GRID_FOLDER {tmpdir / 'missing_grid'}",
                ]
            ),
            encoding="utf-8",
        )

        with self.assertRaises(FileNotFoundError):
            main(["conformer-docking", str(config)])

    def test_conformer_docking_accepts_full_rotation_flag(self):
        parser = build_parser()

        args = parser.parse_args(["conformer-docking", "--full-rotation", "testgrid.in"])

        self.assertTrue(args.full_rotation)

    def test_decompose_accepts_cpp_style_arguments(self):
        tmpdir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(tmpdir, ignore_errors=True))
        output = tmpdir / "annotated.sdf"
        fragments = tmpdir / "fragments.sdf"

        self.assertEqual(
            main(
                [
                    "decompose",
                    "--ligand",
                    str(DATA / "G39.mol2"),
                    "--fragment",
                    str(fragments),
                    "--output",
                    str(output),
                ]
            ),
            0,
        )

        self.assertTrue(output.exists())
        self.assertTrue(fragments.exists())
        self.assertGreater(fragments.read_text(encoding="utf-8").count("$$$$"), 0)
        self.assertIn(">  <fragment_info>", output.read_text(encoding="utf-8"))

    def test_conformer_docking_writes_cpp_style_score_csv(self):
        tmpdir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(tmpdir, ignore_errors=True))
        config = tmpdir / "testgrid.in"
        grid_dir = tmpdir / "grid"
        output = tmpdir / "docked.sdf"
        config.write_text(
            "\n".join(
                [
                    "INNERBOX 2, 2, 2",
                    "OUTERBOX 2, 2, 2",
                    "BOX_CENTER 0.3826, 81.705, 109.195",
                    "SEARCH_PITCH 1, 1, 1",
                    "SCORING_PITCH 1, 1, 1",
                    "MEMORY_SIZE 10",
                    f"RECEPTOR {DATA / '2HU4_A_r.pdbqt'}",
                    f"LIGAND {DATA / 'G39.mol2'}",
                    f"OUTPUT {output}",
                    f"GRID_FOLDER {grid_dir}",
                    "NO_LOCAL_OPT true",
                    "POSES_PER_LIG 1",
                    "POSES_PER_LIG_BEFORE_OPT 1",
                ]
            ),
            encoding="utf-8",
        )
        self.assertEqual(main(["atomgrid-gen", str(config)]), 0)
        self.assertEqual(main(["conformer-docking", str(config)]), 0)

        csvs = list(tmpdir.glob("docked.sdffraggrid__*.csv"))
        self.assertEqual(len(csvs), 1)
        self.assertIn(",", csvs[0].read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
