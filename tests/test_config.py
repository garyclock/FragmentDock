from pathlib import Path
import unittest

from restretto.config import ReuseStrategy, parse_in_file


class ConfigTests(unittest.TestCase):
    def test_parse_reference_config_resolves_main_fields(self):
        conf = parse_in_file(Path("references/restretto/testdata/testgrid.in"))

        self.assertEqual(conf.grid.inner_width.as_tuple(), (10.0, 10.0, 10.0))
        self.assertEqual(conf.grid.outer_width.as_tuple(), (20.0, 20.0, 20.0))
        self.assertEqual(conf.grid.center.as_tuple(), (0.3826, 81.705, 109.195))
        self.assertEqual(conf.grid.search_pitch.as_tuple(), (1.0, 1.0, 1.0))
        self.assertEqual(conf.grid.score_pitch.as_tuple(), (0.25, 0.25, 0.25))
        self.assertEqual(conf.mem_size, 8000)
        self.assertEqual(conf.receptor_file, "testdata/2HU4_A_r.pdbqt")
        self.assertEqual(conf.ligand_files, ["testdata/G39.mol2", "testdata/conformers.mol2"])
        self.assertEqual(conf.output_file, "testdata/G39_docked.sdf")
        self.assertEqual(conf.grid_folder, "testdata/grid")
        self.assertIs(conf.reuse_grid, ReuseStrategy.OFFLINE)
        self.assertFalse(conf.no_local_opt)


    def test_parse_config_rejects_score_only_and_local_only(self):
        config = Path("test_bad_config.in")
        self.addCleanup(lambda: config.exists() and config.unlink())
        config.write_text(
            "\n".join(
                [
                    "OUTERBOX 2, 2, 2",
                    "INNERBOX 1, 1, 1",
                    "BOX_CENTER 0, 0, 0",
                    "SEARCH_PITCH 1, 1, 1",
                    "SCORING_PITCH 0.5, 0.5, 0.5",
                    "MEMORY_SIZE 1",
                    "SCORE_ONLY true",
                    "LOCAL_ONLY true",
                ]
            ),
            encoding="utf-8",
        )

        with self.assertRaisesRegex(ValueError, "SCORE_ONLY.*LOCAL_ONLY"):
            parse_in_file(config)


if __name__ == "__main__":
    unittest.main()
