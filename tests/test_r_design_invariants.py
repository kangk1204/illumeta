"""Executable checks for R design-matrix invariants."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import textwrap
import unittest


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYZE_R = os.path.join(BASE_DIR, "r_scripts", "analyze.R")


def extract_r_function(source: str, name: str) -> str:
    marker = f"{name} <- function"
    start = source.index(marker)
    brace = source.index("{", start)
    depth = 0
    for idx in range(brace, len(source)):
        char = source[idx]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return source[start : idx + 1]
    raise ValueError(f"Could not extract R function: {name}")


class RDesignInvariantTests(unittest.TestCase):
    def run_r(self, code: str) -> None:
        if shutil.which("Rscript") is None:
            self.skipTest("Rscript is not installed")
        with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
            handle.write(code)
            script_path = handle.name
        try:
            result = subprocess.run(
                ["Rscript", script_path],
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
        finally:
            os.unlink(script_path)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)

    def test_rank_pruning_preserves_group_contrast_columns(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        drop_fn = extract_r_function(source, "drop_linear_dependencies")

        self.run_r(
            textwrap.dedent(
                f"""
                {drop_fn}
                mat <- cbind(
                  Control = c(1, 1, 0, 0),
                  Case = c(0, 0, 1, 1),
                  ConfoundedControl = c(1, 1, 0, 0),
                  NumericCovariate = c(1, 2, 3, 4)
                )
                res <- drop_linear_dependencies(mat, c("Control", "Case"))
                if (!all(c("Control", "Case") %in% colnames(res$mat))) {{
                  stop("group contrast columns were not preserved")
                }}
                if ("ConfoundedControl" %in% colnames(res$mat)) {{
                  stop("confounded covariate was retained instead of the group column")
                }}
                if (!("ConfoundedControl" %in% res$dropped)) {{
                  stop("confounded covariate was not reported as dropped")
                }}
                """
            )
        )

    def test_covariate_collinearity_tolerates_all_na_correlations(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        map_fn = extract_r_function(source, "map_mm_columns_to_vars")
        collinearity_fn = extract_r_function(source, "compute_covariate_collinearity")

        self.run_r(
            textwrap.dedent(
                f"""
                {map_fn}
                {collinearity_fn}
                targets <- data.frame(
                  Constant = c(1, 1, 1, 1),
                  Alternating = c(0, 1, 0, 1),
                  NumericCovariate = c(1, 2, 3, 4)
                )
                res <- compute_covariate_collinearity(
                  targets,
                  c("Constant", "Alternating", "NumericCovariate")
                )
                if (nrow(res) != 3) stop("unexpected row count")
                constant_row <- res[res$Variable == "Constant", ]
                if (nrow(constant_row) != 1 || !is.na(constant_row$Max_Correlation)) {{
                  stop("all-NA correlations should be represented as NA")
                }}
                """
            )
        )


if __name__ == "__main__":
    unittest.main()
