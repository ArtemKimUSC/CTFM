# CT-FM custom annotation-panel helpers

This folder contains helper scripts for building a custom CT-FM annotation panel
and the adjusted LD-score correlation matrix used as `R.txt`.

The workflow mirrors the original project notes:

1. Convert BED regulatory regions into thin S-LDSC `.annot.gz` files.
2. Run `ldsc.py --l2 --thin-annot` to compute annotation LD scores.
3. Collect target annotation LD scores into matrix A.
4. Collect Whole/background/baseline LD scores into matrix B.
5. Regress A on B, scale residualized LD scores, and compute annotation correlations.

## Scripts

- `scripts/create_thin_annot_from_bed.py`: replacement for the original
  `create_annot.pl`; creates one thin `.annot` file from BED intervals and a
  chromosome-specific PLINK `.bim`.
- `scripts/create_ldscores_from_bed.sh`: loops over chromosomes, creates
  `.annot.gz` files, and calls LDSC to produce `.l2.ldscore.gz` files.
- `scripts/make_whole_bed.py`: merges all target BED intervals into one Whole
  background BED.
- `scripts/collect_ldscore_matrices.py`: concatenates target and background
  LD-score files into per-chromosome A and B matrices, checking CHR/SNP/BP row
  alignment.
- `scripts/make_adjusted_ldscore_correlation.R`: residualizes A on B and writes
  the adjusted annotation correlation matrix.

## Minimal command order

```bash
python ctfm_custom_annotation_protocol/scripts/make_whole_bed.py \
  --bed-list data/custom_annotations/custom_beds.txt \
  --out data/custom_annotations/beds/Whole.bed

bash ctfm_custom_annotation_protocol/scripts/create_ldscores_from_bed.sh \
  --bed-list data/custom_annotations/custom_beds.txt \
  --out-dir data/custom_annotations/ldscores/custom_panel \
  --name custom_panel \
  --ref-bfile-prefix data/1000G_EUR_Phase3/plink_files/1000G.EUR.QC \
  --ldsc-py /path/to/ldsc.py \
  --print-snps data/1000G_EUR_Phase3/w_hm3.txt

bash ctfm_custom_annotation_protocol/scripts/create_ldscores_from_bed.sh \
  --bed data/custom_annotations/beds/Whole.bed \
  --out-dir data/custom_annotations/ldscores/whole \
  --name Whole \
  --ref-bfile-prefix data/1000G_EUR_Phase3/plink_files/1000G.EUR.QC \
  --ldsc-py /path/to/ldsc.py \
  --print-snps data/1000G_EUR_Phase3/w_hm3.txt

python ctfm_custom_annotation_protocol/scripts/collect_ldscore_matrices.py \
  --annotation-list data/custom_annotations/target_ldscore_items.txt \
  --background-list data/custom_annotations/background_ldscore_items.txt \
  --out-dir data/custom_annotations/matrices \
  --chromosomes 1-22

Rscript ctfm_custom_annotation_protocol/scripts/make_adjusted_ldscore_correlation.R \
  --a-dir data/custom_annotations/matrices/A \
  --b-dir data/custom_annotations/matrices/B \
  --out data/custom_annotations/R_custom.txt \
  --chromosomes 1-22
```

## Manifest examples

`custom_beds.txt`:

```text
/path/to/Osteoblast.bed    Osteoblast
/path/to/CD4_T_cells.bed   CD4_T_cells
/path/to/Monocytes.bed     Monocytes
```

`target_ldscore_items.txt`:

```text
data/custom_annotations/ldscores/custom_panel/custom_panel
```

`background_ldscore_items.txt`:

```text
data/custom_annotations/ldscores/whole/Whole
data/1000G_EUR_Phase3/baseline_v1.2/baseline
```

All LD-score files must be computed using the same PLINK reference panel,
genome build, and LDSC `--print-snps` file.

The matrix collection script requires Python with `pandas`. The correlation
script requires R with `data.table`.
