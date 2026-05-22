"""Tests for smoke-pipeline job parsing."""

from pathlib import Path

from scripts.run_smoke_pipeline import load_jobs, safe_job_name


def test_safe_job_name_blocks_path_traversal():
    assert safe_job_name("../../outside") == "outside"
    assert safe_job_name("nested/path\\name") == "nested_path_name"
    assert safe_job_name("  ") == "job"


def test_load_jobs_sanitizes_derived_name(tmp_path):
    jobs_path = tmp_path / "jobs.tsv"
    jobs_path.write_text(
        "config\tgroup_con\tgroup_test\tname\n"
        "configs/example.tsv\tControl\tCase\t../../outside\n",
        encoding="utf-8",
    )

    jobs = load_jobs(Path(jobs_path))

    assert jobs[0]["raw_name"] == "../../outside"
    assert jobs[0]["name"] == "outside"
