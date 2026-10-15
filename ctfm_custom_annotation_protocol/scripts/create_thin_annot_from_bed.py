#!/usr/bin/env python3
"""Create a thin S-LDSC annotation file from BED interval files and a PLINK BIM.

The output matches the "thin" annotation format used by LDSC:
one header row with annotation names, followed by one 0/1 row per SNP in the BIM
file. Coordinates are interpreted as standard BED coordinates: 0-based,
half-open intervals in the same genome build as the BIM file.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


@dataclass(frozen=True)
class BedAnnotation:
    path: Path
    name: str


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom


def sanitize_name(name: str) -> str:
    for suffix in [".bed.gz", ".bed", ".narrowPeak.gz", ".narrowPeak", ".broadPeak.gz", ".broadPeak", ".gz"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    name = re.sub(r"[^A-Za-z0-9_.+-]+", "_", name)
    name = name.strip("_")
    if not name:
        raise ValueError("Could not derive a non-empty annotation name")
    return name


def looks_like_path(token: str) -> bool:
    return (
        "/" in token
        or token.endswith((".bed", ".bed.gz", ".narrowPeak", ".narrowPeak.gz", ".broadPeak", ".broadPeak.gz"))
    )


def parse_manifest_line(line: str) -> BedAnnotation | None:
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    fields = line.split()
    if len(fields) == 1:
        bed_path = Path(fields[0])
        name = sanitize_name(bed_path.name)
    else:
        first, second = fields[0], fields[1]
        if looks_like_path(first) and not looks_like_path(second):
            bed_path = Path(first)
            name = sanitize_name(second)
        elif looks_like_path(second) and not looks_like_path(first):
            bed_path = Path(second)
            name = sanitize_name(first)
        else:
            bed_path = Path(first)
            name = sanitize_name(second)
    return BedAnnotation(path=bed_path, name=name)


def read_bed_manifest(path: Path) -> list[BedAnnotation]:
    annotations: list[BedAnnotation] = []
    with open_text(path) as handle:
        for raw in handle:
            entry = parse_manifest_line(raw)
            if entry is not None:
                annotations.append(entry)
    return annotations


def collect_annotations(args: argparse.Namespace) -> list[BedAnnotation]:
    annotations: list[BedAnnotation] = []
    if args.bed_list:
        annotations.extend(read_bed_manifest(Path(args.bed_list)))
    for bed in args.bed or []:
        bed_path = Path(bed)
        annotations.append(BedAnnotation(path=bed_path, name=sanitize_name(bed_path.name)))
    if not annotations:
        raise SystemExit("ERROR: provide at least one --bed or --bed-list")

    seen: set[str] = set()
    for annot in annotations:
        if annot.name in seen:
            raise SystemExit(f"ERROR: duplicate annotation name after sanitization: {annot.name}")
        seen.add(annot.name)
    return annotations


def read_bim_positions(bim: Path) -> tuple[str, list[int]]:
    chrom: str | None = None
    positions: list[int] = []
    with bim.open("r") as handle:
        for raw in handle:
            fields = raw.split()
            if len(fields) < 4:
                continue
            this_chrom = normalize_chrom(fields[0])
            if chrom is None:
                chrom = this_chrom
            elif chrom != this_chrom:
                raise SystemExit(
                    f"ERROR: BIM file {bim} contains multiple chromosomes ({chrom}, {this_chrom}). "
                    "Run this script separately for each chromosome."
                )
            positions.append(int(fields[3]))
    if chrom is None or not positions:
        raise SystemExit(f"ERROR: no SNP positions found in {bim}")
    return chrom, positions


def read_intervals_for_chrom(bed: Path, chrom: str) -> list[tuple[int, int]]:
    intervals: list[tuple[int, int]] = []
    with open_text(bed) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith(("#", "track", "browser")):
                continue
            fields = raw.split()
            if len(fields) < 3:
                continue
            if normalize_chrom(fields[0]) != chrom:
                continue
            start = int(fields[1])
            end = int(fields[2])
            if end > start:
                intervals.append((start, end))
    intervals.sort()
    return intervals


def mark_overlaps(positions: Iterable[int], intervals: list[tuple[int, int]]) -> list[int]:
    calls: list[int] = []
    idx = 0
    n_intervals = len(intervals)
    for pos in positions:
        while idx < n_intervals and intervals[idx][1] < pos:
            idx += 1
        if idx < n_intervals and intervals[idx][0] < pos <= intervals[idx][1]:
            calls.append(1)
        else:
            calls.append(0)
    return calls


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a thin S-LDSC .annot file from BED intervals and one chromosome BIM file."
    )
    parser.add_argument("--bim", required=True, help="PLINK BIM file for one chromosome")
    parser.add_argument("--out", required=True, help="Output thin .annot file")
    parser.add_argument("--bed-list", help="Manifest with one BED path per row, optionally with an annotation name")
    parser.add_argument("--bed", action="append", help="BED file; may be provided more than once")
    args = parser.parse_args()

    annotations = collect_annotations(args)
    chrom, positions = read_bim_positions(Path(args.bim))

    columns: list[list[int]] = []
    for annot in annotations:
        intervals = read_intervals_for_chrom(annot.path, chrom)
        calls = mark_overlaps(positions, intervals)
        columns.append(calls)
        print(
            f"{annot.name}: {sum(calls)} / {len(calls)} SNPs overlap intervals on chromosome {chrom}",
            file=sys.stderr,
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        handle.write("\t".join(annot.name for annot in annotations) + "\n")
        for row_idx in range(len(positions)):
            handle.write("\t".join(str(col[row_idx]) for col in columns) + "\n")


if __name__ == "__main__":
    main()
