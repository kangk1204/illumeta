"""Tests for --group-map sample exclusion.

Public GEO series routinely carry more arms than a given contrast uses: an MCI
arm alongside AD and control, a second assay subset, a second tissue. Before
this feature the only way to run a two-arm contrast over such a series was to
hand-edit configure.tsv, because auto-group required every row to resolve to the
control or test label and refused to leave any row empty. That manual step was
never recorded anywhere, which is precisely what makes a published sample count
impossible to reproduce from the deposited metadata.

These tests pin the exclusion path: the right samples are dropped, the drop is
reported rather than silent, and the guards that made auto-group trustworthy in
the first place still fire.
"""

from __future__ import annotations

import csv
import sys
import unittest
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

BASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE_DIR))

import illumeta  # noqa: E402


HEADERS = ["SampleID", "Basename", "primary_group", "disease state:ch1"]

# Three arms of two samples each, mirroring the shape of a real AD/MCI/control series.
ARMS = {
    "cognitively healthy status;cognitively normal": "CTL",
    "late status of cognitive decline;severe dementia": "CASE",
    "early status of cognitive decline;early phase": "MCI",
}


def _write_config(path: Path, arms=ARMS, per_arm: int = 2) -> Path:
    rows = []
    n = 0
    for label in arms:
        for _ in range(per_arm):
            n += 1
            rows.append([f"GSM{n:07d}", f"basename_{n}", "", label])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(HEADERS)
        writer.writerows(rows)
    return path


def _counts(path: str) -> dict[str, int]:
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return dict(Counter(row["primary_group"] for row in rows))


class GroupMapParsingTests(unittest.TestCase):
    def test_exclusion_tokens_map_to_the_sentinel(self):
        for token in ("exclude", "drop", "omit", "skip", "none", "-", "NA", "Excluded"):
            mapping = illumeta.parse_group_map(f"arm_x={token}", "Control", "AD")
            self.assertEqual(
                mapping[illumeta.normalize_group_value("arm_x")],
                illumeta.GROUP_MAP_EXCLUDE,
                msg=f"token {token!r} should mark the arm for exclusion",
            )

    def test_normal_targets_are_unaffected(self):
        mapping = illumeta.parse_group_map("a=Control,b=AD", "Control", "AD")
        self.assertEqual(mapping[illumeta.normalize_group_value("a")], "Control")
        self.assertEqual(mapping[illumeta.normalize_group_value("b")], "AD")

    def test_apply_group_mapping_propagates_the_sentinel(self):
        mapping = illumeta.parse_group_map("arm_x=exclude", "Control", "AD")
        self.assertEqual(
            illumeta.apply_group_mapping("arm_x", mapping, "Control", "AD"),
            illumeta.GROUP_MAP_EXCLUDE,
        )


class AutoGroupExclusionTests(unittest.TestCase):
    def test_third_arm_is_dropped_and_reported(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            cfg = _write_config(tmp / "configure.tsv")
            out = tmp / "configure_autogroup.tsv"
            path, info = illumeta.auto_group_config(
                str(cfg), group_con="CTL", group_test="CASE",
                group_map=(
                    "cognitively healthy status/cognitively normal=CTL,"
                    "late status of cognitive decline/severe dementia=CASE,"
                    "early status of cognitive decline/early phase=exclude"
                ),
                output_path=str(out), overwrite=True,
            )
            self.assertEqual(_counts(path), {"CTL": 2, "CASE": 2})
            self.assertEqual(info["excluded_rows"], 2)
            self.assertEqual(info["total_rows"], 4)
            # The drop must be attributable to a specific arm, not just counted.
            self.assertEqual(
                list(info["excluded_labels"]),
                ["early status of cognitive decline;early phase"],
            )
            # The excluded samples must not survive into the written config.
            with open(path, encoding="utf-8", newline="") as handle:
                text = handle.read()
            self.assertNotIn("early phase", text)
            self.assertNotIn(illumeta.GROUP_MAP_EXCLUDE, text)

    def test_without_exclusion_the_third_arm_still_raises(self):
        """The fail-closed guard must survive the feature: an unmapped arm is an error."""
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            cfg = _write_config(tmp / "configure.tsv")
            with self.assertRaises(ValueError) as ctx:
                illumeta.auto_group_config(
                    str(cfg), group_con="CTL", group_test="CASE",
                    output_path=str(tmp / "out.tsv"), overwrite=True,
                )
            self.assertIn("outside control/test groups", str(ctx.exception))

    def test_excluding_every_arm_is_an_error(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            cfg = _write_config(tmp / "configure.tsv")
            with self.assertRaises(ValueError) as ctx:
                illumeta.auto_group_config(
                    str(cfg), group_con="CTL", group_test="CASE",
                    group_map=(
                        "cognitively healthy status/cognitively normal=exclude,"
                        "late status of cognitive decline/severe dementia=exclude,"
                        "early status of cognitive decline/early phase=exclude"
                    ),
                    output_path=str(tmp / "out.tsv"), overwrite=True,
                )
            self.assertIn("excluded every sample", str(ctx.exception))

    def test_no_exclusion_leaves_counts_and_metadata_untouched(self):
        """A two-arm series must behave exactly as before the feature."""
        two_arm = {
            "cognitively healthy status;cognitively normal": "CTL",
            "late status of cognitive decline;severe dementia": "CASE",
        }
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            cfg = _write_config(tmp / "configure.tsv", arms=two_arm, per_arm=3)
            path, info = illumeta.auto_group_config(
                str(cfg), group_con="CTL", group_test="CASE",
                group_map=(
                    "cognitively healthy status/cognitively normal=CTL,"
                    "late status of cognitive decline/severe dementia=CASE"
                ),
                output_path=str(tmp / "out.tsv"), overwrite=True,
            )
            self.assertEqual(_counts(path), {"CTL": 3, "CASE": 3})
            self.assertEqual(info["excluded_rows"], 0)
            self.assertEqual(info["excluded_labels"], {})


if __name__ == "__main__":
    unittest.main()
