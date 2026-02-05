import os
import tempfile
import unittest

from illumeta import ensure_cross_reactive_lists


class CrossReactiveListTests(unittest.TestCase):
    def _write(self, path, content):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(content)

    def test_build_from_local_sources(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sources_dir = os.path.join(tmpdir, "_sources")
            self._write(
                os.path.join(sources_dir, "48639-non-specific-probes-Illumina450k.csv"),
                "TargetID\ncg00000029\ncg00000108\nBADID\n",
            )
            self._write(
                os.path.join(sources_dir, "HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt"),
                "cg00000236\nch.1.1\njunk\n",
            )
            self._write(
                os.path.join(sources_dir, "13059_2016_1066_MOESM1_ESM.csv"),
                "X\ncg00000029\ncg00000108\n",
            )
            self._write(
                os.path.join(sources_dir, "1-s2.0-S221359601630071X-mmc2.txt"),
                "CpG\ncg00000236\n",
            )
            self._write(
                os.path.join(sources_dir, "1-s2.0-S221359601630071X-mmc3.txt"),
                "CpG\ncg00000363\n",
            )

            created, status = ensure_cross_reactive_lists(tmpdir, allow_network=False)
            self.assertTrue(created)
            self.assertEqual(status, "built_local")

            out_450k = os.path.join(tmpdir, "cross_reactive_450K.tsv")
            out_epic = os.path.join(tmpdir, "cross_reactive_EPIC.tsv")
            self.assertTrue(os.path.exists(out_450k))
            self.assertTrue(os.path.exists(out_epic))

            with open(out_450k, "r", encoding="utf-8") as handle:
                lines = [line.strip() for line in handle if line.strip()]
            self.assertEqual(lines[0], "CpG")
            self.assertIn("cg00000029", lines)
            self.assertIn("cg00000108", lines)
            self.assertIn("cg00000236", lines)
            self.assertIn("ch.1.1", lines)
            self.assertNotIn("BADID", lines)

    def test_existing_outputs_skip(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            lines_450k = ["CpG"] + [f"cg00000{i:03d}" for i in range(10)]
            lines_epic = ["CpG"] + [f"cg00001{i:03d}" for i in range(10)]
            self._write(os.path.join(tmpdir, "cross_reactive_450K.tsv"), "\n".join(lines_450k) + "\n")
            self._write(os.path.join(tmpdir, "cross_reactive_EPIC.tsv"), "\n".join(lines_epic) + "\n")
            self._write(os.path.join(tmpdir, "cross_reactive_EPICv2.tsv"), "\n".join(lines_epic) + "\n")
            created, status = ensure_cross_reactive_lists(tmpdir, allow_network=False)
            self.assertFalse(created)
            self.assertEqual(status, "exists")

    def test_missing_sources_no_network(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            created, status = ensure_cross_reactive_lists(tmpdir, allow_network=False)
            self.assertFalse(created)
            self.assertEqual(status, "missing_sources")


if __name__ == "__main__":
    unittest.main()
