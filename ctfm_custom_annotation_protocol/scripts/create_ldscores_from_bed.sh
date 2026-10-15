#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Create S-LDSC thin annotations and LD scores from BED intervals.

Required:
  --out-dir DIR                 Output directory for .annot.gz and .l2.ldscore.gz files
  --name NAME                   Panel name; used as the output prefix when --bed is provided
  --ref-bfile-prefix PREFIX     PLINK reference prefix, e.g. /path/1000G.EUR.QC
  --ldsc-py PATH                Path to ldsc.py
  --print-snps PATH             SNP list used by LDSC --print-snps, e.g. w_hm3.txt

Provide either:
  --bed BED                     One BED file; may be repeated
  --bed-list FILE               Manifest with BED paths, optionally with annotation names.
                                Each manifest row is written as a separate LDSC prefix.

Optional:
  --create-annot-py PATH        Path to create_thin_annot_from_bed.py
  --chromosomes "1 2 ... 22"    Chromosomes to process (default: 1..22)
  --ld-wind-cm FLOAT            LDSC LD window in cM (default: 1)

Example:
  bash create_ldscores_from_bed.sh \
    --bed-list custom_beds.txt \
    --out-dir data/custom_annotations/ldscores/custom_panel \
    --name custom_panel \
    --ref-bfile-prefix data/1000G_EUR_Phase3/plink_files/1000G.EUR.QC \
    --ldsc-py /path/to/ldsc.py \
    --print-snps data/1000G_EUR_Phase3/w_hm3.txt

  This creates one prefix per annotation, for example:
    data/custom_annotations/ldscores/custom_panel/CD19_B.1.annot.gz
    data/custom_annotations/ldscores/custom_panel/CD19_B.1.l2.ldscore.gz
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CREATE_ANNOT_PY="$SCRIPT_DIR/create_thin_annot_from_bed.py"
CHROMOSOMES="$(seq 1 22)"
LD_WIND_CM="1"
BED_ARGS=()
BED_LIST=""
OUT_DIR=""
NAME=""
REF_BFILE_PREFIX=""
LDSC_PY=""
PRINT_SNPS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bed)
      BED_ARGS+=(--bed "$2")
      shift 2
      ;;
    --bed-list)
      BED_LIST="$2"
      shift 2
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --name)
      NAME="$2"
      shift 2
      ;;
    --ref-bfile-prefix)
      REF_BFILE_PREFIX="$2"
      shift 2
      ;;
    --ldsc-py)
      LDSC_PY="$2"
      shift 2
      ;;
    --print-snps)
      PRINT_SNPS="$2"
      shift 2
      ;;
    --create-annot-py)
      CREATE_ANNOT_PY="$2"
      shift 2
      ;;
    --chromosomes)
      CHROMOSOMES="$2"
      shift 2
      ;;
    --ld-wind-cm)
      LD_WIND_CM="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$OUT_DIR" || -z "$NAME" || -z "$REF_BFILE_PREFIX" || -z "$LDSC_PY" || -z "$PRINT_SNPS" ]]; then
  usage >&2
  exit 2
fi

if [[ -z "$BED_LIST" && ${#BED_ARGS[@]} -eq 0 ]]; then
  echo "ERROR: provide --bed or --bed-list" >&2
  exit 2
fi

mkdir -p "$OUT_DIR"

run_ldsc_for_prefix() {
  local annot_args=("$@")
  local prefix="${annot_args[0]}"
  annot_args=("${annot_args[@]:1}")

  for CHR in $CHROMOSOMES; do
    local annot="$OUT_DIR/$prefix.$CHR.annot"
    echo "[create_ldscores_from_bed] $prefix chromosome $CHR"

    python "$CREATE_ANNOT_PY" \
      --bim "$REF_BFILE_PREFIX.$CHR.bim" \
      --out "$annot" \
      "${annot_args[@]}"

    gzip -f "$annot"

    python "$LDSC_PY" \
      --l2 \
      --bfile "$REF_BFILE_PREFIX.$CHR" \
      --ld-wind-cm "$LD_WIND_CM" \
      --annot "$annot.gz" \
      --thin-annot \
      --out "$OUT_DIR/$prefix.$CHR" \
      --print-snps "$PRINT_SNPS" \
      > "$OUT_DIR/$prefix.$CHR.ldsc.log" 2>&1
  done
}

if [[ -n "$BED_LIST" ]]; then
  manifest_tmp="$(mktemp "${TMPDIR:-/tmp}/ctfm_bed_manifest.XXXXXX")"
  python - "$BED_LIST" > "$manifest_tmp" <<'PY'
import gzip
import re
import sys
from pathlib import Path

def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

def sanitize_name(name):
    for suffix in [".bed.gz", ".bed", ".narrowPeak.gz", ".narrowPeak", ".broadPeak.gz", ".broadPeak", ".gz"]:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    name = re.sub(r"[^A-Za-z0-9_.+-]+", "_", name).strip("_")
    if not name:
        raise SystemExit("ERROR: could not derive a non-empty annotation name")
    return name

def looks_like_path(token):
    return "/" in token or token.endswith((".bed", ".bed.gz", ".narrowPeak", ".narrowPeak.gz", ".broadPeak", ".broadPeak.gz"))

seen = set()
with open_text(Path(sys.argv[1])) as handle:
    for raw in handle:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) == 1:
            bed = fields[0]
            name = sanitize_name(Path(bed).name)
        else:
            first, second = fields[0], fields[1]
            if looks_like_path(first) and not looks_like_path(second):
                bed = first
                name = sanitize_name(second)
            elif looks_like_path(second) and not looks_like_path(first):
                bed = second
                name = sanitize_name(first)
            else:
                bed = first
                name = sanitize_name(second)
        if name in seen:
            raise SystemExit(f"ERROR: duplicate annotation name after sanitization: {name}")
        seen.add(name)
        print(f"{name}\t{bed}")
PY

  cut -f1 "$manifest_tmp" > "$OUT_DIR/annotation_names.txt"

  while IFS=$'\t' read -r annot_name bed_path; do
    run_ldsc_for_prefix "$annot_name" --bed "$bed_path"
  done < "$manifest_tmp"

  rm -f "$manifest_tmp"
else
  run_ldsc_for_prefix "$NAME" "${BED_ARGS[@]}"
fi

