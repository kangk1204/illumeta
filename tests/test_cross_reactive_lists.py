import os
import tempfile
import unittest

from illumeta import ensure_cross_reactive_lists, resolve_cross_reactive_dir


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


    def test_resolve_cross_reactive_dir_relative_config_yaml(self):
        """Relative --config_yaml should resolve relative to configure.tsv dir, not CWD."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create project structure:
            #   tmpdir/project/configure.tsv
            #   tmpdir/project/my_config.yaml  (has local_dir: custom_lists)
            #   tmpdir/project/custom_lists/    (target directory)
            project_dir = os.path.join(tmpdir, "project")
            os.makedirs(project_dir, exist_ok=True)

            config_path = os.path.join(project_dir, "configure.tsv")
            with open(config_path, "w") as f:
                f.write("Basename\tprimary_group\n")

            yaml_path = os.path.join(project_dir, "my_config.yaml")
            with open(yaml_path, "w") as f:
                f.write("cross_reactive:\n  local_dir: custom_lists\n")

            custom_dir = os.path.join(project_dir, "custom_lists")
            os.makedirs(custom_dir, exist_ok=True)

            # Resolve with a RELATIVE config_yaml path
            result = resolve_cross_reactive_dir(config_path, config_yaml_path="my_config.yaml")

            # Should resolve to project_dir/custom_lists, NOT CWD/custom_lists
            self.assertEqual(os.path.realpath(result), os.path.realpath(custom_dir))

    def test_resolve_cross_reactive_dir_absolute_config_yaml(self):
        """Absolute --config_yaml should work as-is."""
        with tempfile.TemporaryDirectory() as tmpdir:
            project_dir = os.path.join(tmpdir, "project")
            os.makedirs(project_dir, exist_ok=True)

            config_path = os.path.join(project_dir, "configure.tsv")
            with open(config_path, "w") as f:
                f.write("Basename\tprimary_group\n")

            yaml_path = os.path.join(project_dir, "my_config.yaml")
            with open(yaml_path, "w") as f:
                f.write("cross_reactive:\n  local_dir: custom_lists\n")

            custom_dir = os.path.join(project_dir, "custom_lists")
            os.makedirs(custom_dir, exist_ok=True)

            # Absolute path should work directly
            result = resolve_cross_reactive_dir(config_path, config_yaml_path=yaml_path)
            self.assertEqual(os.path.realpath(result), os.path.realpath(custom_dir))


if __name__ == "__main__":
    unittest.main()
