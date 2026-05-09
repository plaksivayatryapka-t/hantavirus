#!/usr/bin/env python3
"""Fetch FASTA records for Nextstrain/Auspice tree tips with NCBI accessions."""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
JSON_FILES = [
    REPO_ROOT / "data/nextstrain/groups/hodcroftlab/andv/S.json",
    REPO_ROOT / "data/nextstrain/groups/hodcroftlab/andv/M.json",
    REPO_ROOT / "data/nextstrain/groups/hodcroftlab/andv/L.json",
]
OUT_DIR = REPO_ROOT / "data/fasta/nextstrain-tips"
NEXTSTRAIN_DATASET_URL = "https://nextstrain.org/groups/hodcroftlab/andv"
NCBI_SOURCE_URL = NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_EMAIL = None
BATCH_SIZE = 100
SLEEP_SECONDS = 0.34
RETRIES = 5
INCLUDE_NONSTANDARD_ACCESSIONS = False

ACCESSION_RE = re.compile(r"^[A-Z]{1,4}_?[0-9]+(?:\.[0-9]+)?$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract leaf accessions from one or more Nextstrain Auspice JSON files "
            "and fetch their nucleotide FASTA records from NCBI E-utilities."
        )
    )
    parser.add_argument(
        "json_files",
        nargs="*",
        type=Path,
        help="Optional Auspice JSON file override(s). Defaults to the JSON_FILES constant.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=OUT_DIR,
        help="Directory for output FASTA and manifest files. Default: %(default)s",
    )
    parser.add_argument(
        "--email",
        default=NCBI_EMAIL,
        help="Email address to pass to NCBI E-utilities.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=BATCH_SIZE,
        help="Number of accessions per NCBI efetch request. Default: %(default)s",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=SLEEP_SECONDS,
        help="Seconds to sleep between NCBI requests. Default: %(default)s",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=RETRIES,
        help="Retries per NCBI request. Default: %(default)s",
    )
    parser.add_argument(
        "--include-nonstandard",
        action="store_true",
        help="Include accessions that do not look like normal INSDC/NCBI accessions.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write only the manifest; do not fetch FASTA records.",
    )
    return parser.parse_args()


def iter_tips(node: dict):
    children = node.get("children") or []
    if not children:
        yield node
        return
    for child in children:
        yield from iter_tips(child)


def node_attr_value(node: dict, key: str):
    value = (node.get("node_attrs") or {}).get(key)
    if isinstance(value, dict):
        return value.get("value")
    return value


def load_tip_accessions(json_file: Path, include_nonstandard: bool) -> list[dict[str, str]]:
    with json_file.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    tree = data.get("tree")
    if not isinstance(tree, dict):
        raise ValueError(f"{json_file} does not contain an Auspice v2 tree object")

    rows = []
    seen = set()
    segment = json_file.stem
    for tip in iter_tips(tree):
        accession = node_attr_value(tip, "accession")
        if not accession or not isinstance(accession, str):
            continue
        accession = accession.strip()
        if not include_nonstandard and not ACCESSION_RE.match(accession):
            continue
        key = (segment, accession)
        if key in seen:
            continue
        seen.add(key)
        rows.append(
            {
                "segment": segment,
                "accession": accession,
                "tip_name": tip.get("name", ""),
                "country": node_attr_value(tip, "country") or "",
            }
        )
    return rows


def chunks(items: list[str], size: int):
    for index in range(0, len(items), size):
        yield items[index : index + size]


def fetch_fasta_batch(
    accessions: list[str],
    email: str | None,
    retries: int,
) -> str:
    form = {
        "db": "nuccore",
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
        "tool": "hantavirus-nextstrain-tip-fetcher",
    }
    if email:
        form["email"] = email
    body = urllib.parse.urlencode(form).encode("utf-8")

    last_error = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(NCBI_EFETCH_URL, data=body, method="POST")
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                text = response.read().decode("utf-8")
            if not text.startswith(">"):
                raise RuntimeError(f"NCBI returned non-FASTA response: {text[:200]!r}")
            return text
        except (urllib.error.URLError, TimeoutError, RuntimeError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(min(2**attempt, 30))
    raise RuntimeError(f"failed to fetch {len(accessions)} accession(s): {last_error}")


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("segment\taccession\ttip_name\tcountry\n")
        for row in rows:
            handle.write(
                f"{row['segment']}\t{row['accession']}\t{row['tip_name']}\t{row['country']}\n"
            )


def write_index(path: Path, json_files: list[Path]) -> None:
    relative_json_files = []
    for json_file in json_files:
        try:
            relative_json_files.append(json_file.relative_to(REPO_ROOT))
        except ValueError:
            relative_json_files.append(json_file)

    lines = [
        f"source: {NEXTSTRAIN_DATASET_URL}",
        f"fasta_source: {NCBI_SOURCE_URL}",
        "inputs:",
    ]
    lines.extend(f"- {json_file}" for json_file in relative_json_files)
    lines.extend(
        [
            "outputs:",
            "- S.fasta",
            "- M.fasta",
            "- L.fasta",
            "- manifest.tsv",
            "",
            "FASTA records are fetched from NCBI nuccore accessions found on leaf nodes in the Nextstrain Auspice JSON files.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    if args.batch_size < 1:
        raise SystemExit("--batch-size must be >= 1")

    all_rows = []
    json_files = args.json_files or JSON_FILES
    include_nonstandard = args.include_nonstandard or INCLUDE_NONSTANDARD_ACCESSIONS
    for json_file in json_files:
        all_rows.extend(load_tip_accessions(json_file, include_nonstandard))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest = args.out_dir / "manifest.tsv"
    write_manifest(manifest, all_rows)
    write_index(args.out_dir / "index.md", json_files)

    by_segment: dict[str, list[str]] = {}
    for row in all_rows:
        by_segment.setdefault(row["segment"], []).append(row["accession"])

    print(f"wrote {len(all_rows)} accession rows to {manifest}", file=sys.stderr)
    if args.dry_run:
        return 0

    for segment, accessions in sorted(by_segment.items()):
        fasta_path = args.out_dir / f"{segment}.fasta"
        with fasta_path.open("w", encoding="utf-8") as handle:
            for batch in chunks(accessions, args.batch_size):
                text = fetch_fasta_batch(batch, args.email, args.retries)
                handle.write(text)
                if not text.endswith("\n"):
                    handle.write("\n")
                time.sleep(args.sleep)
        print(f"wrote {len(accessions)} {segment} FASTA record(s) to {fasta_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
