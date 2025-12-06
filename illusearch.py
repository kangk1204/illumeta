#!/usr/bin/env python3
"""
Search GEO (GSE) for Illumina methylation datasets (450k/850k/950k) that mention IDAT files
and save results to a TSV. Uses NCBI E-utilities (esearch/esummary); pass --email to follow
NCBI etiquette. Optionally checks the GEO FTP supplement listing for actual .idat/.idat.gz/RAW.tar files.
"""

import argparse
import csv
import re
import sys
import time
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

try:
    import requests
except ImportError:
    print("Error: Python package 'requests' is missing.")
    print("Install with: python3 -m pip install requests")
    sys.exit(1)

GEO_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PLATFORM_ACCESSIONS = ["GPL13534", "GPL21145", "GPL34372"]  # 450k, EPIC(850k), EPIC v2 (950k)
# Match .idat/.idat.gz anywhere in the listing, or RAW.tar bundles (no end-of-string anchoring).
IDAT_REGEX = re.compile(r"(\.idat(\.gz)?|RAW\.tar)", re.IGNORECASE)
TOOL_NAME = "illusearch-idat"
ESUMMARY_BATCH_SIZE = 150
PLATFORM_LABELS = {"GPL13534": "450k", "GPL21145": "850k", "GPL34372": "950k"}
SUPPL_RETRY = 3
SUPPL_RETRY_DELAY = 0.5


def build_search_term(keywords: str) -> str:
    if not keywords.strip():
        raise ValueError("keywords must be a non-empty string")
    platform_query = " OR ".join([f"{p}[Accession]" for p in PLATFORM_ACCESSIONS])
    keyword_clause = f"({keywords})"
    term = f"{keyword_clause} AND (idat[All Fields]) AND ({platform_query}) AND gse[Entry Type]"
    return term


def eutils_request(path: str, params: Dict[str, str], email: Optional[str], sleep_s: float) -> requests.Response:
    if email:
        params["email"] = email
    params["tool"] = TOOL_NAME
    url = f"{GEO_EUTILS}/{path}"
    resp = requests.get(url, params=params, timeout=20)
    time.sleep(sleep_s)  # be polite to NCBI
    resp.raise_for_status()
    return resp


def search_gse_ids(term: str, email: Optional[str], sleep_s: float, retmax: int) -> List[str]:
    print(f"[*] Searching GEO GSE with term: {term}")
    params = {
        "db": "gds",
        "term": term,
        "retmax": str(retmax),
        "retmode": "xml",
    }
    resp = eutils_request("esearch.fcgi", params, email, sleep_s)
    root = ET.fromstring(resp.text)
    ids = [elem.text for elem in root.findall(".//Id") if elem.text]
    print(f"[*] Found {len(ids)} candidate GSE records")
    return ids


def fetch_summaries(ids: List[str], email: Optional[str], sleep_s: float) -> List[Dict[str, str]]:
    if not ids:
        return []
    print(f"[*] Fetching summaries for {len(ids)} IDs (batch size: {ESUMMARY_BATCH_SIZE})...")
    out: List[Dict[str, str]] = []
    for start in range(0, len(ids), ESUMMARY_BATCH_SIZE):
        chunk = ids[start : start + ESUMMARY_BATCH_SIZE]
        id_str = ",".join(chunk)
        params = {"db": "gds", "id": id_str, "retmode": "xml"}
        resp = eutils_request("esummary.fcgi", params, email, sleep_s)
        root = ET.fromstring(resp.text)
        before = len(out)
        for doc in root.findall(".//DocSum"):
            rec: Dict[str, str] = {}
            for item in doc.findall("Item"):
                name = item.get("Name")
                if name:
                    if len(item):
                        children = [(child.text or "").strip() for child in item.findall("Item")]
                        value = ";".join([c for c in children if c])
                    else:
                        value = (item.text or "").strip()
                    rec[name.lower()] = value

            accession = rec.get("accession")
            if accession and accession.startswith("GSE"):
                # Normalize platform ID(s) and infer array generation
                raw_gpl = rec.get("gpl", "")
                tokens = [t for t in re.split(r"[;,\s]+", raw_gpl) if t]
                gpls = []
                for tok in tokens:
                    tok = tok.upper()
                    if tok.startswith("GPL"):
                        gpls.append(tok)
                    elif tok.isdigit():
                        gpls.append(f"GPL{tok}")
                    else:
                        gpls.append(tok)
                platform = ";".join(gpls)
                platform_types = sorted({PLATFORM_LABELS[g] for g in gpls if g in PLATFORM_LABELS})
                platform_label = ";".join(platform_types)
                pubmed_ids = re.sub(r"\s+", " ", rec.get("pubmedids", "")).strip()

                out.append(
                    {
                        "gse_id": accession,
                        "title": rec.get("title", ""),
                        "platform": platform,
                        "platform_type": platform_label,
                        "sample_count": rec.get("n_samples", ""),
                        "pubmed_id": pubmed_ids,
                    }
                )
        parsed = len(out) - before
        print(f"    - Parsed {parsed} from batch {start // ESUMMARY_BATCH_SIZE + 1}")
    print(f"[*] Parsed {len(out)} GSE summaries total")
    return out


def geo_suppl_url(gse_id: str) -> str:
    prefix = gse_id[:-3] + "nnn"  # GSE12345 -> GSE12nnn
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_id}/suppl/"


def fetch_with_retry(url: str, method: str = "get") -> Optional[requests.Response]:
    last_resp: Optional[requests.Response] = None
    for attempt in range(1, SUPPL_RETRY + 1):
        try:
            resp = requests.request(method=method, url=url, timeout=15, allow_redirects=True)
            last_resp = resp
            if resp.status_code == 200:
                return resp
        except requests.RequestException:
            last_resp = None
        if attempt < SUPPL_RETRY:
            time.sleep(SUPPL_RETRY_DELAY)
    return last_resp


def has_idat_in_suppl(gse_id: str) -> str:
    url = geo_suppl_url(gse_id)
    resp = fetch_with_retry(url, method="get")
    if resp is None or resp.status_code != 200:
        return "error"
    if IDAT_REGEX.search(resp.text):
        return "yes"

    # Fallback: check for bundled RAW tar directly (HEAD to reduce bandwidth)
    raw_url = f"{url}{gse_id}_RAW.tar"
    raw_resp = fetch_with_retry(raw_url, method="head")
    if raw_resp is not None and raw_resp.status_code == 200:
        return "yes"
    return "no"


def enrich_with_suppl_check(rows: List[Dict[str, str]], check_suppl: bool, sleep_s: float) -> None:
    if not check_suppl:
        return
    print("[*] Checking GEO supplementary directories for actual IDAT/RAW files...")
    for idx, row in enumerate(rows, start=1):
        gse_id = row["gse_id"]
        sys.stdout.write(f"    [{idx}/{len(rows)}] {gse_id}: ")
        sys.stdout.flush()
        found = has_idat_in_suppl(gse_id)
        row["suppl_has_idat"] = found
        if found == "yes":
            print("YES")
        elif found == "no":
            print("no")
        else:
            print("error")
        time.sleep(sleep_s)


def write_tsv(rows: List[Dict[str, str]], path: str) -> None:
    if not rows:
        print("[!] No results to write.")
        return
    fieldnames = ["gse_id", "title", "platform", "platform_type", "sample_count", "pubmed_id", "suppl_has_idat"]
    for row in rows:
        row.setdefault("suppl_has_idat", "unchecked")
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[*] Saved {len(rows)} records to {path}")


def main():
    parser = argparse.ArgumentParser(
        description="Search GEO for IDAT-backed Illumina methylation GSEs (450k/850k/950k) and export to TSV."
    )
    parser.add_argument(
        "-k",
        "--keywords",
        required=True,
        help="Required keywords to AND with the IDAT/platform filter (e.g., 'breast cancer', 'blood cohort').",
    )
    parser.add_argument("-o", "--output", default="geo_idat_methylation.tsv", help="Output TSV path.")
    parser.add_argument("--email", default=None, help="Your email (NCBI etiquette; simple alternative to an API key).")
    parser.add_argument("--retmax", type=int, default=500, help="Maximum records to fetch (default: 500).")
    parser.add_argument(
        "--check-suppl",
        dest="check_suppl",
        action="store_true",
        default=True,
        help="Fetch GEO FTP supplementary listings to confirm presence of .idat/.idat.gz/RAW.tar. (default: on)",
    )
    parser.add_argument(
        "--no-check-suppl",
        dest="check_suppl",
        action="store_false",
        help="Disable supplementary listing check to speed up.",
    )
    parser.add_argument("--sleep", type=float, default=0.5, help="Sleep seconds between E-utility requests (0.5 is gentle without an API key).")
    args = parser.parse_args()

    if args.email:
        print(f"[*] Using contact email for NCBI requests: {args.email}")

    keywords = args.keywords.strip()
    if not keywords:
        parser.error("--keywords is required and cannot be empty.")

    term = build_search_term(keywords)
    ids = search_gse_ids(term, args.email, args.sleep, args.retmax)
    summaries = fetch_summaries(ids, args.email, args.sleep)
    enrich_with_suppl_check(summaries, args.check_suppl, args.sleep)
    write_tsv(summaries, args.output)
    print("[*] Done.")


if __name__ == "__main__":
    main()
