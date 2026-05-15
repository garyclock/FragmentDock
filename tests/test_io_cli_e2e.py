from pathlib import Path
import contextlib
import io
import shutil
import tempfile
import unittest

from restretto.cli import main
from restretto.io import CHEM_BACKEND, read_molecules


DATA = Path("references/restretto/testdata")


class IoCliE2ETests(unittest.TestCase):
    def test_read_reference_mol2_and_pdbqt(self):
        self.assertEqual(CHEM_BACKEND, "openbabel")
        ligands = read_molecules(DATA / "G39.mol2")
        receptor = read_molecules(DATA / "2HU4_A_r.pdbqt")

        self.assertEqual(len(ligands), 1)
        self.assertEqual(ligands[0].size, 23)
        self.assertGreaterEqual(len(ligands[0].bonds), 23)
        self.assertGreater(receptor[0].size, 100)
        self.assertGreater(len(receptor[0].bonds), 100)


    def test_atomgrid_score_decompose_and_docking_cli(self):
        tmpdir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(tmpdir, ignore_errors=True))
        config = tmpdir / "testgrid.in"
        grid_dir = tmpdir / "grid"
        output = tmpdir / "docked.sdf"
        config.write_text(
            "\n".join(
                [
                    "INNERBOX 4, 4, 4",
                    "OUTERBOX 6, 6, 6",
                    "BOX_CENTER 0.3826, 81.705, 109.195",
                    "SEARCH_PITCH 1, 1, 1",
                    "SCORING_PITCH 1, 1, 1",
                    "MEMORY_SIZE 10",
                    f"RECEPTOR {DATA / '2HU4_A_r.pdbqt'}",
                    f"LIGAND {DATA / 'G39.mol2'}",
                    f"OUTPUT {output}",
                    f"GRID_FOLDER {grid_dir}",
                    "NO_LOCAL_OPT true",
                    "POSES_PER_LIG 2",
                ]
            ),
            encoding="utf-8",
        )

        self.assertEqual(main(["atomgrid-gen", str(config)]), 0)
        self.assertTrue(any(grid_dir.glob("*.grid")))

        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(["score-only", str(config)]), 0)
            self.assertEqual(main(["intraenergy-only", str(config)]), 0)
        self.assertEqual(main(["decompose", str(config), "--output", str(tmpdir / "fragments.sdf")]), 0)
        self.assertEqual(main(["conformer-docking", str(config)]), 0)
        self.assertTrue(output.exists())
        self.assertGreaterEqual(output.read_text(encoding="utf-8").count("$$$$"), 1)


if __name__ == "__main__":
    unittest.main()
