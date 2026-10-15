#!/usr/bin/env python3
"""Merge a BED manifest into one union ("Whole") BED annotation."""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from typing import TextIO


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def normalize_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    c = normalize_chrom(chrom)
    if c.isdigit():
        return (int(c), "")
    return (10_000, c)


def looks_like_path(token: str) -> bool:
    return "/" in token or token.endswith((".bed", ".bed.gz", ".narrowPeak", ".narrowPeak.gz", ".broadPeak", ".broadPeak.gz"))


def read_manifest(path: Path) -> list[Path]:
    beds: list[Path] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) == 1:
                beds.append(Path(fields[0]))
            elif looks_like_path(fields[0]):
                beds.append(Path(fields[0]))
            else:
                beds.append(Path(fields[1]))
    return beds


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a merged union BED from a BED manifest.")
    parser.add_argument("--bed-list", required=True, help="Manifest with BED paths")
    parser.add_argument("--out", required=True, help="Output merged BED file")
    parser.add_argument(
        "--chr-prefix",
        choices=["preserve", "add", "remove"],
        default="add",
        help="How chromosome names should be written in the output BED",
    )
    args = parser.parse_args()

    intervals_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    beds = read_manifest(Path(args.bed_list))
    if not beds:
        raise SystemExit("ERROR: no BED files found in manifest")

    for bed in beds:
        with open_text(bed) as handle:
            for raw in handle:
                if not raw.strip() or raw.startswith(("#", "track", "browser")):
                    continue
                fields = raw.split()
                if len(fields) < 3:
                    continue
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                if end <= start:
                    continue
                intervals_by_chrom[normalize_chrom(chrom)].append((start, end))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    n_out = 0
    with out.open("w") as handle:
        for chrom in sorted(intervals_by_chrom, key=chrom_sort_key):
            intervals = sorted(intervals_by_chrom[chrom])
            merged: list[tuple[int, int]] = []
            for start, end in intervals:
                if not merged or start > merged[-1][1]:
                    merged.append((start, end))
                else:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            if args.chr_prefix == "add":
                chrom_out = f"chr{chrom}"
            elif args.chr_prefix == "remove":
                chrom_out = chrom
            else:
                chrom_out = chrom if chrom.lower().startswith("chr") else f"chr{chrom}"
            for start, end in merged:
                handle.write(f"{chrom_out}\t{start}\t{end}\n")
                n_out += 1

    print(f"Wrote {n_out} merged intervals to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
