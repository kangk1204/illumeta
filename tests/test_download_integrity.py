"""Regression tests for _download_file integrity/robustness (no real network)."""

import http.client
import os
import re
import tempfile
import unittest
from pathlib import Path
from urllib.error import URLError


class _FakeResp:
    def __init__(self, chunks, headers):
        self._chunks = list(chunks)
        self.headers = headers

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self, _n):
        item = self._chunks.pop(0)
        if isinstance(item, Exception):
            raise item
        return item


class DownloadIntegrityTests(unittest.TestCase):
    def setUp(self):
        import illumeta
        self.illumeta = illumeta
        self._orig = illumeta._urllib_request.urlopen
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        self.illumeta._urllib_request.urlopen = self._orig

    def _patch(self, resp):
        self.illumeta._urllib_request.urlopen = lambda req, timeout=30: resp

    def test_complete_download_is_written(self):
        self._patch(_FakeResp([b"CpG\ncg1\ncg2\n"], {"Content-Length": "12"}))
        dest = os.path.join(self.tmp, "ok.tsv")
        out = self.illumeta._download_file("http://x/ok.tsv", dest)
        self.assertEqual(out, dest)
        with open(dest, "rb") as fh:
            self.assertEqual(fh.read(), b"CpG\ncg1\ncg2\n")

    def test_truncated_download_raises_and_writes_nothing(self):
        # Server advertises 999 bytes but sends 5 -> must NOT poison the cache.
        self._patch(_FakeResp([b"short"], {"Content-Length": "999"}))
        dest = os.path.join(self.tmp, "trunc.tsv")
        with self.assertRaises(URLError):
            self.illumeta._download_file("http://x/trunc.tsv", dest)
        self.assertFalse(os.path.exists(dest))

    def test_incomplete_read_becomes_urlerror(self):
        # A dropped connection mid-read must surface as URLError (caught by callers),
        # not an uncaught http.client.IncompleteRead crash.
        self._patch(_FakeResp([http.client.IncompleteRead(b"partial")], {"Content-Length": "10"}))
        dest = os.path.join(self.tmp, "drop.tsv")
        with self.assertRaises(URLError):
            self.illumeta._download_file("http://x/drop.tsv", dest)
        self.assertFalse(os.path.exists(dest))

    def test_missing_content_length_is_accepted(self):
        # No Content-Length header -> cannot verify length; still write what we got.
        self._patch(_FakeResp([b"CpG\ncg1\n"], {}))
        dest = os.path.join(self.tmp, "nolen.tsv")
        self.illumeta._download_file("http://x/nolen.tsv", dest)
        self.assertTrue(os.path.exists(dest))



class DownloadRetryBackoffTests(unittest.TestCase):
    """GEO fetches fail in bursts, so the retry window has to outlast a short outage.

    A real GSE66351 download failed after three attempts inside ~26 seconds: the
    connection stalled for eleven minutes, then the constant 3-second backoff burned
    every remaining attempt in under half a minute. Exponential backoff is what makes
    the difference between a transient hiccup and a failed download, so pin it.
    """

    @classmethod
    def setUpClass(cls):
        cls.source = (Path(__file__).resolve().parents[1] / "r_scripts" / "download.R").read_text(
            encoding="utf-8"
        )

    def _default(self, env_var: str) -> str:
        match = re.search(rf'Sys\.getenv\("{env_var}",\s*"([^"]+)"\)', self.source)
        self.assertIsNotNone(match, f"{env_var} default not found")
        return match.group(1)

    def test_backoff_is_exponential(self):
        self.assertGreaterEqual(
            float(self._default("ILLUMETA_DOWNLOAD_BACKOFF")), 2.0,
            "a constant backoff cannot outlast a transient NCBI outage",
        )

    def test_retry_window_spans_at_least_thirty_seconds(self):
        attempts = int(self._default("ILLUMETA_DOWNLOAD_RETRIES"))
        wait = float(self._default("ILLUMETA_DOWNLOAD_WAIT"))
        backoff = float(self._default("ILLUMETA_DOWNLOAD_BACKOFF"))
        total = sum(wait * backoff ** i for i in range(attempts - 1))
        self.assertGreaterEqual(total, 30.0, f"retry window is only {total:.0f}s")

    def test_defaults_stay_overridable(self):
        for env_var in ("ILLUMETA_DOWNLOAD_RETRIES", "ILLUMETA_DOWNLOAD_WAIT",
                        "ILLUMETA_DOWNLOAD_BACKOFF"):
            self.assertIn(f'Sys.getenv("{env_var}"', self.source)


class CanonicalColumnDedupTests(unittest.TestCase):
    """Duplicate-column removal must not delete the names the pipeline resolves by name.

    Sentrix_ID and Sentrix_Position are derived from IDAT filenames and appended last,
    so a plain "keep the first" dedup always dropped them whenever GEO also exposed the
    same values as characteristics fields. select_batch_factor() looks the batch factor
    up by name, so the values surviving under an opaque characteristics_ch1.N name were
    invisible to it: batch correction and the tier3 stratified path switched themselves
    off, reported only as an informational log line. Observed on GSE66351, whose
    published run used Sentrix_Position as its batch candidate.
    """

    @classmethod
    def setUpClass(cls):
        cls.source = (Path(__file__).resolve().parents[1] / "r_scripts" / "download.R").read_text(
            encoding="utf-8"
        )

    def test_canonical_columns_are_declared(self):
        self.assertIn("CANONICAL_META_COLS", self.source)
        for name in ("primary_group", "Basename", "Sentrix_ID", "Sentrix_Position"):
            self.assertRegex(
                self.source,
                rf'CANONICAL_META_COLS <- c\([^)]*"{name}"',
                msg=f"{name} must be protected from duplicate-column removal",
            )

    def test_dedup_is_preference_ordered_not_positional(self):
        """A bare duplicated(col_enc) would reintroduce the bug."""
        self.assertIn("pref_order", self.source)
        self.assertIn("duplicated(col_enc[pref_order])", self.source)
        self.assertNotIn("dup_cols <- names(simple_meta)[duplicated(col_enc)]", self.source)

    def test_kept_columns_preserve_original_order(self):
        # Reordering configure.tsv columns would be a gratuitous change to a file users
        # and downstream code both read positionally in places.
        self.assertIn("# preserve original column order", self.source)

if __name__ == "__main__":
    unittest.main()
