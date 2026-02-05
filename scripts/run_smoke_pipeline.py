#!/usr/bin/env python3
import argparse
import csv
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path


REQUIRED_FIELDS = ("config", "group_con", "group_test")


def load_jobs(path: Path):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Job file is empty or missing headers.")
        for field in REQUIRED_FIELDS:
            if field not in reader.fieldnames:
                raise ValueError(f"Missing required column: {field}")
        jobs = []
        for row in reader:
            job = {k: (v or "").strip() for k, v in row.items()}
            if not job["config"]:
                continue
            job["config"] = str(Path(job["config"]).expanduser())
            job["name"] = job.get("name") or Path(job["config"]).stem
            jobs.append(job)
    if not jobs:
        raise ValueError("No valid jobs found.")
    return jobs


def run_job(job, args, out_dir: Path, log_dir: Path):
    output = job.get("output") or str(out_dir / job["name"])
    log_path = log_dir / f"{job['name']}.log"

    cmd = [
        args.python,
        args.illumeta,
        "analysis",
        "--config",
        job["config"],
        "--group_con",
        job["group_con"],
        "--group_test",
        job["group_test"],
        "-o",
        output,
    ]
    if args.extra_args:
        cmd.extend(shlex.split(args.extra_args))
    if job.get("extra_args"):
        cmd.extend(shlex.split(job["extra_args"]))

    env = os.environ.copy()
    if args.r_libs_user:
        env["R_LIBS_USER"] = args.r_libs_user
    env["ILLUMETA_TRACEBACK"] = "1"

    start = time.time()
    with log_path.open("w") as logf:
        proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, env=env)
    duration = time.time() - start

    summary_path = Path(output) / "summary.json"
    dashboard_path = Path(f"{output}_index.html")
    return {
        "name": job["name"],
        "config": job["config"],
        "group_con": job["group_con"],
        "group_test": job["group_test"],
        "output": output,
        "exit_code": str(proc.returncode),
        "duration_sec": f"{duration:.1f}",
        "summary_exists": str(summary_path.is_file()),
        "dashboard_exists": str(dashboard_path.is_file()),
        "log_path": str(log_path),
    }


def write_report(rows, report_path: Path):
    fields = [
        "name",
        "config",
        "group_con",
        "group_test",
        "output",
        "exit_code",
        "duration_sec",
        "summary_exists",
        "dashboard_exists",
        "log_path",
    ]
    with report_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def main():
    parser = argparse.ArgumentParser(
        description="Run IlluMeta analysis smoke tests from a TSV job file."
    )
    parser.add_argument("--jobs", required=True, help="TSV with columns: config, group_con, group_test, [name], [output], [extra_args].")
    parser.add_argument("--out-dir", default="benchmarks/smoke_runs", help="Output root for smoke runs.")
    parser.add_argument("--illumeta", default="illumeta.py", help="Path to illumeta.py.")
    parser.add_argument("--python", default=sys.executable, help="Python executable to run illumeta.")
    parser.add_argument("--r-libs-user", default="", help="R_LIBS_USER to inject into runs.")
    parser.add_argument("--extra-args", default="", help="Extra args appended to every run.")
    parser.add_argument("--stop-on-failure", action="store_true", help="Stop after first failed run.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    log_dir = out_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    jobs = load_jobs(Path(args.jobs))
    rows = []
    for job in jobs:
        row = run_job(job, args, out_dir, log_dir)
        rows.append(row)
        if args.stop_on_failure and row["exit_code"] != "0":
            break

    report_path = out_dir / "smoke_report.tsv"
    write_report(rows, report_path)


if __name__ == "__main__":
    main()
