#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Run S-LDSC with a custom CT-FM annotation panel.

Required:
  --trait TRAIT                 Trait basename matching sumstats_ready/TRAIT.sumstats.gz
  --custom-dir DIR              Directory containing annots.txt and R.txt, e.g. data/SLDSC_custom

Provide either:
  --ldcts FILE                  Pre-built LDSC .ldcts file with one row per annotation

or all of:
  --target-prefix-template STR  Template for one-target LD-score prefixes; use {ANNOT}
  --whole-prefix PREFIX         Whole/background LD-score prefix

Optional:
  --out-root DIR                S-LDSC output root (default: out/SLDSC_custom)
  --sumstats-dir DIR            Summary-statistic directory (default: sumstats_ready)
  --ref-dir DIR                 1000G reference directory (default: data/1000G_EUR_Phase3)
  --ldsc-py PATH_OR_CMD         LDSC command or ldsc.py path (default: ldsc.py)
  --block-size INT              Number of annotations per S-LDSC block (default: 200)
  --mode local|slurm|print      Execute locally, submit with sbatch, or print commands (default: local)
  --mem MEM                     Slurm memory when --mode slurm (default: 40G)
  --time TIME                   Slurm time when --mode slurm (default: 0-24:00:00)

Examples:
  bash launch_sldsc_custom.sh \
    --trait YOUR_TRAIT \
    --custom-dir data/SLDSC_custom \
    --ldcts data/SLDSC_custom/custom.ldcts \
    --ldsc-py /project2/gazal_569/DATA/ldsc/ldsc-2.0.1/ldsc.py

  bash launch_sldsc_custom.sh \
    --trait YOUR_TRAIT \
    --custom-dir data/SLDSC_custom \
    --target-prefix-template 'data/custom_annotations/ldscores/{ANNOT}/{ANNOT}.' \
    --whole-prefix data/custom_annotations/ldscores/whole/Whole. \
    --ldsc-py /project2/gazal_569/DATA/ldsc/ldsc-2.0.1/ldsc.py
EOF
}

TRAIT=""
CUSTOM_DIR=""
LDCTS=""
TARGET_PREFIX_TEMPLATE=""
WHOLE_PREFIX=""
OUT_ROOT="out/SLDSC_custom"
SUMSTATS_DIR="sumstats_ready"
REF_DIR="data/1000G_EUR_Phase3"
LDSC_PY="ldsc.py"
BLOCK_SIZE="200"
MODE="local"
SLURM_MEM="40G"
SLURM_TIME="0-24:00:00"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --trait) TRAIT="$2"; shift 2 ;;
    --custom-dir) CUSTOM_DIR="$2"; shift 2 ;;
    --ldcts) LDCTS="$2"; shift 2 ;;
    --target-prefix-template) TARGET_PREFIX_TEMPLATE="$2"; shift 2 ;;
    --whole-prefix) WHOLE_PREFIX="$2"; shift 2 ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --sumstats-dir) SUMSTATS_DIR="$2"; shift 2 ;;
    --ref-dir) REF_DIR="$2"; shift 2 ;;
    --ldsc-py) LDSC_PY="$2"; shift 2 ;;
    --block-size) BLOCK_SIZE="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --mem) SLURM_MEM="$2"; shift 2 ;;
    --time) SLURM_TIME="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$TRAIT" || -z "$CUSTOM_DIR" ]]; then
  usage >&2
  exit 2
fi

if [[ -z "$LDCTS" && ( -z "$TARGET_PREFIX_TEMPLATE" || -z "$WHOLE_PREFIX" ) ]]; then
  echo "ERROR: provide --ldcts, or both --target-prefix-template and --whole-prefix" >&2
  exit 2
fi

if [[ "$MODE" != "local" && "$MODE" != "slurm" && "$MODE" != "print" ]]; then
  echo "ERROR: --mode must be local, slurm, or print" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

CUSTOM_DIR="${CUSTOM_DIR%/}"
OUTDIR="$OUT_ROOT/$TRAIT"
BLOCK_DIR="$CUSTOM_DIR/sldsc_blocks"
SUMSTAT="$SUMSTATS_DIR/$TRAIT.sumstats.gz"
BASELINE="$REF_DIR/baseline_v1.2/baseline."
WEIGHTS="$REF_DIR/weights/weights.hm3_noMHC."

mkdir -p "$OUTDIR" "$BLOCK_DIR"

if [[ ! -f "$SUMSTAT" ]]; then
  echo "ERROR: missing summary statistics: $SUMSTAT" >&2
  exit 1
fi

if [[ ! -f "$CUSTOM_DIR/annots.txt" ]]; then
  echo "ERROR: missing custom annotation list: $CUSTOM_DIR/annots.txt" >&2
  exit 1
fi

if [[ -z "$LDCTS" ]]; then
  LDCTS="$CUSTOM_DIR/custom.ldcts"
  python - "$CUSTOM_DIR/annots.txt" "$TARGET_PREFIX_TEMPLATE" "$WHOLE_PREFIX" "$LDCTS" <<'PY'
from pathlib import Path
import sys

annots_path, target_template, whole_prefix, out_path = sys.argv[1:5]

names = []
with open(annots_path) as handle:
    for raw in handle:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        first = line.split()[0]
        if first == "Name":
            continue
        names.append(first)

if not names:
    raise SystemExit(f"ERROR: no annotations found in {annots_path}")

out = Path(out_path)
out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w") as handle:
    for name in names:
        target_prefix = target_template.replace("{ANNOT}", name)
        handle.write(f"{name}\t{target_prefix},{whole_prefix}\n")

print(f"Wrote {len(names)} rows to {out}")
PY
fi

if [[ ! -f "$LDCTS" ]]; then
  echo "ERROR: missing custom .ldcts file: $LDCTS" >&2
  exit 1
fi

python - "$LDCTS" "$BLOCK_DIR" "$BLOCK_SIZE" <<'PY'
from pathlib import Path
import sys

ldcts, block_dir, block_size = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3])
rows = []
with open(ldcts) as handle:
    for raw in handle:
        line = raw.strip()
        if line and not line.startswith("#"):
            rows.append(line)

if not rows:
    raise SystemExit(f"ERROR: no rows found in {ldcts}")

block_dir.mkdir(parents=True, exist_ok=True)
for old in block_dir.glob("sldsc_custom_*"):
    if old.is_file():
        old.unlink()

for block_idx, start in enumerate(range(0, len(rows), block_size)):
    block = rows[start:start + block_size]
    path = block_dir / f"sldsc_custom_{block_idx:02d}"
    path.write_text("\n".join(block) + "\n")
    print(path)
PY

run_one_block() {
  local block_file="$1"
  local block_id="$2"
  local outfile="$OUTDIR/$TRAIT.$block_id"
  local tmp_cts="$outfile.cts"

  # Convert repository-relative data/ paths to absolute paths for LDSC.
  sed "s@data/@$REPO_ROOT/data/@g" "$block_file" > "$tmp_cts"

  local ldsc_cmd
  if [[ -f "$LDSC_PY" ]]; then
    ldsc_cmd="python \"$LDSC_PY\""
  else
    ldsc_cmd="$LDSC_PY"
  fi

  local cmd="$ldsc_cmd --h2-cts \"$SUMSTAT\" --ref-ld-chr \"$BASELINE\" --ref-ld-chr-cts \"$tmp_cts\" --w-ld-chr \"$WEIGHTS\" --out \"$outfile\""

  if [[ "$MODE" == "print" ]]; then
    echo "$cmd"
  elif [[ "$MODE" == "slurm" ]]; then
    sbatch --mem="$SLURM_MEM" -t "$SLURM_TIME" --wrap="$cmd; rm -f \"$tmp_cts\""
  else
    eval "$cmd"
    rm -f "$tmp_cts"
  fi
}

block_id=0
for block_file in "$BLOCK_DIR"/sldsc_custom_*; do
  printf -v block_label "%02d" "$block_id"
  run_one_block "$block_file" "$block_label"
  block_id=$((block_id + 1))
done

echo "Custom S-LDSC launch complete for $TRAIT with $block_id block(s)."

