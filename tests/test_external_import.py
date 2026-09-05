"""Contract tests for `illumeta.py import-external`.

The conversion sits upstream of arithmetic that is already validated against metafor,
so its failure mode is not a wrong number -- it is a right number computed from the
wrong column. These tests therefore concentrate on the ways a table can be misread
silently: a column matched by accident, an SE that was never reported, effects on two
different scales pooled as though they were one.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from illumeta_external import (  # noqa: E402
    ExternalCohortSpec,
    convert_cohort,
    load_external_manifest,
    resolve_columns,
    run_import_external,
)


def write_table(path: Path, header: list[str], rows: list[list[object]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    return path


def champ_like(path: Path, n: int = 5) -> Path:
    """ChAMP returns a limma table: an effect and a t, but no explicit SE column."""
    return write_table(
        path,
        ["CpG", "logFC", "t", "P.Value", "adj.P.Val", "gene", "CHR", "MAPINFO"],
        [[f"cg{i:08d}", 0.10 + i / 100, 2.0 + i, 0.01 * (i + 1), 0.05, "GENE%d" % i, "1", 1000 + i]
         for i in range(n)],
    )


def meffil_like(path: Path, n: int = 5) -> Path:
    """meffil reports the moderated SE directly, so no reconstruction is needed."""
    return write_table(
        path,
        ["CpG", "coefficient", "coefficient.se", "t.statistic", "p.value"],
        [[f"cg{i:08d}", 0.10 + i / 100, 0.05, 2.0 + i, 0.01 * (i + 1)] for i in range(n)],
    )


def read_written(out_dir: Path, filename: str = "Minfi_DMPs_full.csv") -> list[dict[str, str]]:
    with (out_dir / filename).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def base_spec(table: Path, **kw) -> ExternalCohortSpec:
    params = dict(cohort_id="COHORT1", table=table, n_con=10, n_test=10, effect_scale="m")
    params.update(kw)
    return ExternalCohortSpec(**params)


class TestColumnResolution:
    def test_limma_style_headers_need_no_mapping(self, tmp_path):
        resolved, _ = resolve_columns(
            ["CpG", "logFC", "SE", "t", "P.Value"], {}, "C")
        assert resolved["CpG"] == "CpG"
        assert resolved["logFC"] == "logFC"
        assert resolved["SE"] == "SE"

    def test_explicit_mapping_wins_over_autodetection(self):
        # Both 'beta' and 'estimate' would auto-match the effect slot; an explicit
        # choice must be honoured rather than resolved by candidate order.
        resolved, _ = resolve_columns(
            ["probe", "beta", "estimate", "stderr", "pval"],
            {"col_cpg": "probe", "col_effect": "estimate", "col_se": "stderr", "col_pvalue": "pval"},
            "C",
        )
        assert resolved["logFC"] == "estimate"

    def test_explicit_mapping_to_a_missing_column_is_an_error(self):
        # Falling back to auto-detection here would let a typo produce a run that
        # succeeds while pooling a column the user never named.
        with pytest.raises(ValueError, match="not"):
            resolve_columns(["CpG", "logFC", "SE", "P.Value"], {"col_effect": "logFCC"}, "C")

    def test_ambiguous_header_is_refused_rather_than_guessed(self):
        with pytest.raises(ValueError, match="normalise"):
            resolve_columns(["CpG", "logFC", "log_FC", "SE", "P.Value"], {}, "C")

    def test_missing_effect_column_is_an_error(self):
        with pytest.raises(ValueError, match="logFC"):
            resolve_columns(["CpG", "P.Value", "SE"], {"col_effect": ""}, "C")

    def test_no_se_and_no_t_is_refused(self):
        # A P-value cannot supply the weight; accepting the table would mean inventing one.
        with pytest.raises(ValueError, match="standard-error"):
            resolve_columns(["CpG", "logFC", "P.Value"], {}, "C")


class TestConversion:
    def test_champ_style_table_reconstructs_se_from_t(self, tmp_path):
        res = convert_cohort(base_spec(champ_like(tmp_path / "champ.csv")), tmp_path / "out")
        assert res.se_source == "derived_from_t"
        rows = read_written(res.out_dir)
        assert len(rows) == 5
        for row in rows:
            # SE = |logFC / t| is exact for any Wald statistic, so the written t must
            # reproduce the input t up to sign.
            assert math.isclose(abs(float(row["logFC"]) / float(row["SE"])),
                                abs(float(row["t"])), rel_tol=1e-12)

    def test_meffil_style_table_uses_the_reported_se(self, tmp_path):
        res = convert_cohort(base_spec(meffil_like(tmp_path / "meffil.csv")), tmp_path / "out")
        assert res.se_source == "reported"
        assert all(math.isclose(float(r["SE"]), 0.05, rel_tol=1e-12) for r in read_written(res.out_dir))

    def test_effect_scale_is_required_and_validated(self, tmp_path):
        with pytest.raises(ValueError, match="effect_scale"):
            convert_cohort(base_spec(champ_like(tmp_path / "c.csv"), effect_scale=""), tmp_path / "o")
        with pytest.raises(ValueError, match="effect_scale"):
            convert_cohort(base_spec(champ_like(tmp_path / "c2.csv"), effect_scale="logit"), tmp_path / "o2")

    def test_beta_scale_is_accepted_but_flagged(self, tmp_path):
        res = convert_cohort(base_spec(champ_like(tmp_path / "c.csv"), effect_scale="beta"), tmp_path / "o")
        assert any("beta scale" in w for w in res.warnings)

    def test_duplicate_cpgs_are_dropped_and_counted(self, tmp_path):
        table = write_table(
            tmp_path / "dup.csv",
            ["CpG", "logFC", "t", "P.Value"],
            [["cg1", 0.1, 2.0, 0.01], ["cg1", 0.2, 3.0, 0.02], ["cg2", 0.1, 2.0, 0.01]],
        )
        res = convert_cohort(base_spec(table), tmp_path / "out")
        assert res.n_dropped_duplicate == 1
        assert res.n_rows_written == 2

    def test_rows_without_a_usable_effect_or_se_are_dropped(self, tmp_path):
        table = write_table(
            tmp_path / "na.csv",
            ["CpG", "logFC", "t", "P.Value"],
            [["cg1", "NA", "NA", 0.5], ["cg2", 0.1, 0.0, 0.5], ["cg3", 0.1, 2.0, 0.01]],
        )
        res = convert_cohort(base_spec(table), tmp_path / "out")
        assert res.n_dropped_unusable == 2
        assert res.n_rows_written == 1

    def test_a_table_with_nothing_usable_fails_loudly(self, tmp_path):
        table = write_table(tmp_path / "empty.csv", ["CpG", "logFC", "t", "P.Value"],
                            [["cg1", "NA", "NA", 0.5]])
        with pytest.raises(ValueError, match="no usable rows"):
            convert_cohort(base_spec(table), tmp_path / "out")

    def test_missing_delta_beta_is_warned_because_a_gate_depends_on_it(self, tmp_path):
        res = convert_cohort(base_spec(champ_like(tmp_path / "c.csv")), tmp_path / "out")
        assert any("Delta_Beta" in w for w in res.warnings)

    def test_summary_records_that_the_estimates_are_not_illumetas(self, tmp_path):
        spec = base_spec(champ_like(tmp_path / "c.csv"), tool="ChAMP", tool_version="2.40.0")
        res = convert_cohort(spec, tmp_path / "out")
        summary = json.loads((res.out_dir / "summary.json").read_text())
        assert summary["n_con"] == 10 and summary["n_test"] == 10
        assert summary["primary_result_mode"] == "external_import"
        prov = summary["external_import"]
        assert prov["tool"] == "ChAMP" and prov["tool_version"] == "2.40.0"
        assert prov["se_source"] == "derived_from_t"
        assert prov["resolved_columns"]["logFC"] == "logFC"

    def test_external_mode_does_not_trip_the_tier3_variance_warning(self, tmp_path):
        # _tier3_variance_warnings fires on any mode starting with 'tier3_'. An import
        # that named itself that way would attach a within-cohort-meta warning to a
        # cohort that never ran one.
        from illumeta_meta import _tier3_variance_warnings
        res = convert_cohort(base_spec(champ_like(tmp_path / "c.csv")), tmp_path / "out")
        summary = json.loads((res.out_dir / "summary.json").read_text())
        assert _tier3_variance_warnings("C", summary) == []

    def test_branch_alias_is_resolved_to_a_known_file(self, tmp_path):
        res = convert_cohort(base_spec(champ_like(tmp_path / "c.csv"), branch="sesame"), tmp_path / "out")
        assert res.branch == "sesame_strict"
        assert res.table_written.name == "Sesame_DMPs_full.csv"

    def test_unknown_branch_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="branch"):
            convert_cohort(base_spec(champ_like(tmp_path / "c.csv"), branch="nope"), tmp_path / "out")


class TestManifest:
    def _manifest(self, tmp_path, rows: list[dict[str, object]]) -> Path:
        path = tmp_path / "manifest.tsv"
        fields = sorted({k for r in rows for k in r})
        with path.open("w", encoding="utf-8", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        return path

    def test_group_sizes_are_required(self, tmp_path):
        m = self._manifest(tmp_path, [{"cohort": "A", "table": str(champ_like(tmp_path / "a.csv")),
                                       "effect_scale": "m", "n_test": 10}])
        with pytest.raises(ValueError, match="n_con"):
            load_external_manifest(m, tmp_path)

    def test_mixing_effect_scales_across_cohorts_is_refused(self, tmp_path):
        m = self._manifest(tmp_path, [
            {"cohort": "A", "table": str(champ_like(tmp_path / "a.csv")),
             "effect_scale": "m", "n_con": 10, "n_test": 10},
            {"cohort": "B", "table": str(champ_like(tmp_path / "b.csv")),
             "effect_scale": "beta", "n_con": 10, "n_test": 10},
        ])
        with pytest.raises(ValueError, match="effect scale"):
            run_import_external(m, tmp_path / "out", tmp_path, log=lambda *_: None)

    def test_excluded_rows_are_skipped(self, tmp_path):
        m = self._manifest(tmp_path, [
            {"cohort": "A", "table": str(champ_like(tmp_path / "a.csv")),
             "effect_scale": "m", "n_con": 10, "n_test": 10, "include": "1"},
            {"cohort": "B", "table": str(champ_like(tmp_path / "b.csv")),
             "effect_scale": "m", "n_con": 10, "n_test": 10, "include": "no"},
        ])
        assert [s.cohort_id for s in load_external_manifest(m, tmp_path)] == ["A"]

    def test_output_is_readable_by_the_meta_layer(self, tmp_path):
        # The whole point of the command: what it writes must load through the same
        # manifest reader the meta CLI uses, with no special-casing.
        from illumeta_meta import load_meta_manifest
        rows = [{"cohort": g, "table": str(champ_like(tmp_path / f"{g}.csv")),
                 "effect_scale": "m", "n_con": 10, "n_test": 12, "tool": "ChAMP"}
                for g in ("GSE1", "GSE2", "GSE3")]
        m = self._manifest(tmp_path, rows)
        _, meta_manifest = run_import_external(m, tmp_path / "out", tmp_path, log=lambda *_: None)
        cohorts = load_meta_manifest(meta_manifest, tmp_path)
        assert [c.cohort_id for c in cohorts] == ["GSE1", "GSE2", "GSE3"]
        assert all(c.n_con == 10 and c.n_test == 12 for c in cohorts)
        assert all(not c.warnings for c in cohorts)

    def test_end_to_end_pooling_matches_a_hand_computed_fixed_effect(self, tmp_path):
        # Three cohorts with identical effects and SEs: the inverse-variance pooled
        # estimate must be that same effect, and its SE must shrink by sqrt(3).
        from illumeta_meta import load_meta_manifest, _read_branch_records, _random_effect_meta_one
        rows = []
        for g in ("A", "B", "C"):
            t = write_table(tmp_path / f"{g}.csv", ["CpG", "logFC", "SE", "P.Value"],
                            [["cg1", 0.2, 0.05, 0.001]])
            rows.append({"cohort": g, "table": str(t), "effect_scale": "m",
                         "n_con": 10, "n_test": 10})
        m = self._manifest(tmp_path, rows)
        _, meta_manifest = run_import_external(m, tmp_path / "out", tmp_path, log=lambda *_: None)
        cohorts = load_meta_manifest(meta_manifest, tmp_path)
        records, _, _, _ = _read_branch_records(cohorts, "minfi", "Minfi_DMPs_full.csv", False)
        rec = records["cg1"]
        res = _random_effect_meta_one(rec["effects"], rec["ses"], [True] * 3)
        assert math.isclose(res["fixed_effect"], 0.2, rel_tol=1e-9)
        assert math.isclose(res["fixed_se"], 0.05 / math.sqrt(3), rel_tol=1e-9)
        assert res["k"] == 3
