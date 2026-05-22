"""CLI failure paths must return non-zero status codes."""

import os
import subprocess
import sys
import unittest


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ILLUMETA = os.path.join(BASE_DIR, "illumeta.py")


class CliExitCodeTests(unittest.TestCase):
    def run_illumeta(self, *args):
        return subprocess.run(
            [sys.executable, ILLUMETA, *args],
            capture_output=True,
            text=True,
            timeout=30,
        )

    def test_search_empty_keywords_is_nonzero(self):
        result = self.run_illumeta("search", "--keywords", "", "--no-check-suppl")

        self.assertEqual(result.returncode, 2)
        self.assertIn("--keywords is required", result.stderr)

    def test_download_missing_gse_is_nonzero(self):
        result = self.run_illumeta("download")

        self.assertEqual(result.returncode, 2)
        self.assertIn("Example: python illumeta.py download", result.stdout)

    def test_missing_subcommand_is_nonzero(self):
        result = self.run_illumeta()

        self.assertEqual(result.returncode, 2)
        self.assertIn("Available commands", result.stdout)


if __name__ == "__main__":
    unittest.main()
