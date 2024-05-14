
TRAIT=$1
LDCTS=$2
OUT_NUM=$3

SCRIPT_PATH=$(dirname $0)
DATA_PATH=$SCRIPT_PATH/../data

echo $DATA_PATH

SLDSC="python ~/bin/ldsc/ldsc.py";
REF_DIR="$DATA_PATH/1000G_EUR_Phase3"
BASELINE="$REF_DIR/baseline_v1.2/baseline.";
WEIGHTS="--w-ld-chr $REF_DIR/weights/weights.hm3_noMHC.";


# modify SUMSTAT and OUTDIR variables as needed
SUMSTAT="sumstats_ready/$TRAIT.sumstats.gz"
OUTDIR="out/SLDSC/$TRAIT/"

OUTFILE="$OUTDIR/$TRAIT.$OUT_NUM"
mkdir -p $OUTDIR

#sed "s@data@${DATA_PATH}@g" $SCRIPT_PATH/../$LDCTS > $OUTFILE.cts

sed "s@data@${DATA_PATH}@g" $LDCTS > $OUTFILE.cts


python ~/bin/ldsc/ldsc.py \
  --h2-cts $SUMSTAT \
  --ref-ld-chr $BASELINE \
  --ref-ld-chr-cts $OUTFILE.cts \
  $WEIGHTS \
  --out $OUTFILE

rm $OUTFILE.cts
