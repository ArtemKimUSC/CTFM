#!/usr/bin/env python3
"""Collect annotation and background LD scores into per-chromosome matrices.

This replaces brittle shell paste commands with explicit key checks on
CHR/SNP/BP. Each input line may be a directory containing
*.CHR.l2.ldscore.gz files, a prefix for PREFIX.CHR.l2.ldscore.gz, or a template
containing "{chr}".
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path


KEY_COLS = ["CHR", "SNP", "BP"]
pd = None


def read_items(path: Path) -> list[str]:
    items: list[str] = []
    with path.open("r") as handle:
        for raw in handle:
            line = raw.strip()
            if line and not line.startswith("#"):
                items.append(line.split()[0])
    if not items:
        raise SystemExit(f"ERROR: no items found in {path}")
    return items


def parse_chromosomes(value: str) -> list[str]:
    if "-" in value and "," not in value and " " not in value:
        start, end = value.split("-", 1)
        return [str(i) for i in range(int(start), int(end) + 1)]
    return [x for x in re.split(r"[\s,]+", value.strip()) if x]


def resolve_l2_file(item: str, chrom: str) -> Path:
    if "{chr}" in item:
        return Path(item.format(chr=chrom))

    path = Path(item)
    if path.is_dir():
        matches = sorted(path.glob(f"*.{chrom}.l2.ldscore.gz"))
        if len(matches) != 1:
            raise SystemExit(f"ERROR: expected one *.{chrom}.l2.ldscore.gz in {path}, found {len(matches)}")
        return matches[0]

    candidate = Path(f"{item}.{chrom}.l2.ldscore.gz")
    if candidate.exists() or not path.exists():
        return candidate
    return path


def basename_for_item(item: str) -> str:
    path = Path(item.replace("{chr}", "CHR"))
    name = path.name
    for suffix in [".l2.ldscore.gz", ".ldscore.gz", ".gz"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def make_unique(names: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for name in names:
        if name not in seen:
            seen[name] = 0
            out.append(name)
        else:
            seen[name] += 1
            out.append(f"{name}.{seen[name]}")
    return out


def strip_ldscore_suffix(name: str) -> str:
    """Remove the terminal L2 added by LDSC to annotation LD-score columns."""
    return name[:-2] if name.endswith("L2") else name


def read_ldscore(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"ERROR: missing LD-score file: {path}")
    compression = "gzip" if str(path).endswith(".gz") else None
    data = pd.read_csv(path, sep=r"\s+", compression=compression)
    missing = [c for c in KEY_COLS if c not in data.columns]
    if missing:
        raise SystemExit(f"ERROR: {path} is missing key columns: {missing}")
    value_cols = [c for c in data.columns if c not in KEY_COLS]
    if not value_cols:
        raise SystemExit(f"ERROR: {path} has no LD-score value columns")
    return data[KEY_COLS + value_cols]


def append_item(
    matrix: pd.DataFrame | None,
    keys: pd.DataFrame | None,
    item: str,
    chrom: str,
    single_column_name: str,
    strip_l2_suffix: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = resolve_l2_file(item, chrom)
    data = read_ldscore(path)
    item_keys = data[KEY_COLS]
    values = data.drop(columns=KEY_COLS).copy()

    if values.shape[1] == 1 and single_column_name == "basename":
        values.columns = [basename_for_item(item)]
    elif strip_l2_suffix:
        values.columns = [strip_ldscore_suffix(str(c)) for c in values.columns]

    if keys is None:
        keys = item_keys
    elif not keys.equals(item_keys):
        raise SystemExit(
            f"ERROR: SNP rows in {path} do not match the first LD-score file for chromosome {chrom}. "
            "Use the same reference panel and --print-snps list for all annotations."
        )

    if matrix is None:
        matrix = values
    else:
        matrix = pd.concat([matrix, values], axis=1)
    matrix.columns = make_unique(list(matrix.columns))
    return matrix, keys


def build_matrix(
    items: list[str],
    chrom: str,
    single_column_name: str,
    strip_l2_suffix: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix: pd.DataFrame | None = None
    keys: pd.DataFrame | None = None
    for item in items:
        matrix, keys = append_item(matrix, keys, item, chrom, single_column_name, strip_l2_suffix)
    assert matrix is not None and keys is not None
    return matrix, keys


def write_matrix(matrix: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(path, sep="\t", index=False, compression="gzip" if str(path).endswith(".gz") else None)


def main() -> None:
    global pd
    parser = argparse.ArgumentParser(description="Build A and B LD-score matrices for a custom CT-FM R matrix.")
    parser.add_argument("--annotation-list", required=True, help="List of target annotation LD-score dirs/prefixes/templates")
    parser.add_argument("--background-list", required=True, help="List of background LD-score dirs/prefixes/templates")
    parser.add_argument("--out-dir", required=True, help="Output directory containing A/, B/, and keys/")
    parser.add_argument("--chromosomes", default="1-22", help="Chromosomes, e.g. 1-22 or '1 2 3'")
    parser.add_argument(
        "--annotation-single-column-name",
        choices=["basename", "ldscore"],
        default="basename",
        help="How to name one-column target LD-score files",
    )
    parser.add_argument(
        "--background-single-column-name",
        choices=["basename", "ldscore"],
        default="ldscore",
        help="How to name one-column background LD-score files",
    )
    parser.add_argument(
        "--keep-ldscore-l2-suffix",
        action="store_true",
        help="Keep the terminal L2 suffix that LDSC appends to annotation LD-score column names",
    )
    args = parser.parse_args()

    try:
        import pandas as pandas_module
    except ImportError as exc:
        raise SystemExit(
            "ERROR: collect_ldscore_matrices.py requires pandas. "
            "Install it in the active environment, for example with `conda install pandas`."
        ) from exc
    pd = pandas_module

    annotation_items = read_items(Path(args.annotation_list))
    background_items = read_items(Path(args.background_list))
    out_dir = Path(args.out_dir)

    annotation_names: list[str] | None = None
    background_names: list[str] | None = None
    for chrom in parse_chromosomes(args.chromosomes):
        print(f"[collect_ldscore_matrices] chromosome {chrom}", file=sys.stderr)
        strip_l2 = not args.keep_ldscore_l2_suffix
        a_matrix, keys = build_matrix(annotation_items, chrom, args.annotation_single_column_name, strip_l2)
        b_matrix, b_keys = build_matrix(background_items, chrom, args.background_single_column_name, strip_l2)
        if not keys.equals(b_keys):
            raise SystemExit(f"ERROR: target and background SNP rows differ for chromosome {chrom}")

        write_matrix(a_matrix, out_dir / "A" / f"L2.{chrom}.txt.gz")
        write_matrix(b_matrix, out_dir / "B" / f"L2.{chrom}.txt.gz")
        write_matrix(keys, out_dir / "keys" / f"L2.{chrom}.keys.txt.gz")

        if annotation_names is None:
            annotation_names = list(a_matrix.columns)
            background_names = list(b_matrix.columns)
        elif annotation_names != list(a_matrix.columns):
            raise SystemExit(f"ERROR: annotation columns differ on chromosome {chrom}")
        elif background_names != list(b_matrix.columns):
            raise SystemExit(f"ERROR: background columns differ on chromosome {chrom}")

    assert annotation_names is not None and background_names is not None
    (out_dir / "annotation_names.txt").write_text("\n".join(annotation_names) + "\n")
    (out_dir / "background_names.txt").write_text("\n".join(background_names) + "\n")


if __name__ == "__main__":
    main()

