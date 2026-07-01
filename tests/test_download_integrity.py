"""Regression tests for _download_file integrity/robustness (no real network)."""

import http.client
import os
import tempfile
import unittest
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


if __name__ == "__main__":
    unittest.main()
