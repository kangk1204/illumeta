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



class SampleFilterTests(unittest.TestCase):
    """--sample-filter restricts to an assay/tissue subset before grouping.

    --group-map can only exclude values of the grouping column, but GEO series
    routinely put the subset in a different field from the diagnosis. GSE66351
    deposits 190 samples: 128 bulk tissue plus 31 sorted neuron and 31 sorted glia in
    characteristics_ch1, with the diagnosis in characteristics_ch1.1. Reproducing the
    bulk-only contrast without this option means hand-editing configure.tsv.
    """

    HEADERS = ["SampleID", "Basename", "primary_group", "tissue_type", "diagnosis"]

    def _config(self, tmp: Path) -> Path:
        rows, n = [], 0
        for tissue, diagnosis, count in (
            ("bulk", "CTRL", 4), ("bulk", "AD", 6),
            ("Neuron", "CTRL", 2), ("Neuron", "AD", 2),
            ("Glia", "CTRL", 2), ("Glia", "AD", 2),
        ):
            for _ in range(count):
                n += 1
                rows.append([f"GSM{n:07d}", f"basename_{n}", "", tissue, diagnosis])
        path = tmp / "configure.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(self.HEADERS)
            writer.writerows(rows)
        return path

    def test_subset_column_differs_from_group_column(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            path, info = illumeta.auto_group_config(
                str(self._config(tmp)), group_con="Control", group_test="AD",
                group_column="diagnosis", group_map="CTRL=Control,AD=AD",
                sample_filters=["tissue_type=bulk"],
                output_path=str(tmp / "out.tsv"), overwrite=True,
            )
            self.assertEqual(_counts(path), {"Control": 4, "AD": 6})
            entry = info["sample_filters"][0]
            self.assertEqual((entry["kept"], entry["dropped"]), (10, 8))
            self.assertEqual(info["rows_before_filter"], 18)

    def test_multiple_values_and_repeated_filters_combine(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            path, _info = illumeta.auto_group_config(
                str(self._config(tmp)), group_con="Control", group_test="AD",
                group_column="diagnosis", group_map="CTRL=Control,AD=AD",
                sample_filters=["tissue_type=Neuron|Glia", "diagnosis=AD"],
                output_path=str(tmp / "out.tsv"), overwrite=True,
            )
            self.assertEqual(_counts(path), {"AD": 4})

    def test_unknown_column_is_an_error(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            with self.assertRaises(ValueError) as ctx:
                illumeta.auto_group_config(
                    str(self._config(tmp)), group_con="Control", group_test="AD",
                    group_column="diagnosis", sample_filters=["tissue=bulk"],
                    output_path=str(tmp / "out.tsv"), overwrite=True,
                )
            # A typo must fail, not silently analyse every sample.
            self.assertIn("not in configure.tsv", str(ctx.exception))
            self.assertIn("tissue_type", str(ctx.exception))

    def test_filter_matching_nothing_is_an_error(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            with self.assertRaises(ValueError) as ctx:
                illumeta.auto_group_config(
                    str(self._config(tmp)), group_con="Control", group_test="AD",
                    group_column="diagnosis", sample_filters=["tissue_type=cortex"],
                    output_path=str(tmp / "out.tsv"), overwrite=True,
                )
            self.assertIn("matched no samples", str(ctx.exception))

    def test_malformed_filter_is_an_error(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            with self.assertRaises(ValueError):
                illumeta.auto_group_config(
                    str(self._config(tmp)), group_con="Control", group_test="AD",
                    group_column="diagnosis", sample_filters=["tissue_type"],
                    output_path=str(tmp / "out.tsv"), overwrite=True,
                )

class SampleFilterRegexTests(unittest.TestCase):
    """--sample-filter COLUMN~REGEX, and reaching columns configure.tsv dropped.

    GSE105109 deposits 384 samples: every donor appears twice, once bisulfite-treated
    and once oxidative-bisulfite-treated, and the two are distinguishable ONLY through
    the sample title ("entorhinal cortex_bs_1" against "entorhinal cortex_oxbs_1").
    Two things therefore have to work: matching a pattern inside a free-text field, and
    reading a field that building configure.tsv discards as degenerate but that the
    configure_original.tsv snapshot still holds.
    """

    HEADERS = ["geo_accession", "Basename", "primary_group", "diagnosis"]
    FULL_HEADERS = HEADERS + ["title"]

    def _configs(self, tmp: Path) -> Path:
        rows, full, n = [], [], 0
        for assay in ("bs", "oxbs"):
            for diagnosis, count in (("Control", 2), ("AD", 3)):
                for _ in range(count):
                    n += 1
                    acc = f"GSM{n:07d}"
                    rows.append([acc, f"basename_{n}", "", diagnosis])
                    full.append([acc, f"basename_{n}", "", diagnosis,
                                 f"entorhinal cortex_{assay}_{n}"])
        cfg = tmp / "configure.tsv"
        for path, headers, data in ((cfg, self.HEADERS, rows),
                                    (tmp / "configure_original.tsv", self.FULL_HEADERS, full)):
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(headers)
                writer.writerows(data)
        return cfg

    def test_regex_reaches_a_column_only_in_the_snapshot(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            path, info = illumeta.auto_group_config(
                str(self._configs(tmp)), group_con="Control", group_test="AD",
                group_column="diagnosis", sample_filters=[r"title~_bs_"],
                output_path=str(tmp / "out.tsv"), overwrite=True,
            )
            # _bs_ must not also match _oxbs_: the character before "bs" differs.
            self.assertEqual(_counts(path), {"Control": 2, "AD": 3})
            entry = info["sample_filters"][0]
            self.assertEqual(entry["mode"], "regex")
            self.assertEqual((entry["kept"], entry["dropped"]), (5, 5))

    def test_regex_alternation_is_not_split_on_pipe(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            path, _info = illumeta.auto_group_config(
                str(self._configs(tmp)), group_con="Control", group_test="AD",
                group_column="diagnosis", sample_filters=[r"title~_bs_|_oxbs_"],
                output_path=str(tmp / "out.tsv"), overwrite=True,
            )
            self.assertEqual(_counts(path), {"Control": 4, "AD": 6})

    def test_exact_mode_still_splits_on_pipe(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            path, info = illumeta.auto_group_config(
                str(self._configs(tmp)), group_con="Control", group_test="AD",
                group_column="diagnosis", sample_filters=["diagnosis=Control|AD"],
                output_path=str(tmp / "out.tsv"), overwrite=True,
            )
            self.assertEqual(info["sample_filters"][0]["mode"], "exact")
            self.assertEqual(_counts(path), {"Control": 4, "AD": 6})

    def test_invalid_regex_is_an_error(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            with self.assertRaises(ValueError) as ctx:
                illumeta.auto_group_config(
                    str(self._configs(tmp)), group_con="Control", group_test="AD",
                    group_column="diagnosis", sample_filters=[r"title~(unclosed"],
                    output_path=str(tmp / "out.tsv"), overwrite=True,
                )
            self.assertIn("invalid regex", str(ctx.exception))

    def test_column_in_neither_file_is_still_an_error(self):
        with TemporaryDirectory() as raw:
            tmp = Path(raw)
            with self.assertRaises(ValueError) as ctx:
                illumeta.auto_group_config(
                    str(self._configs(tmp)), group_con="Control", group_test="AD",
                    group_column="diagnosis", sample_filters=[r"nosuchfield~x"],
                    output_path=str(tmp / "out.tsv"), overwrite=True,
                )
            self.assertIn("not in configure.tsv", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
