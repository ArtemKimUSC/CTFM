#eval "$(conda shell.bash hook)"
#module purge
#module load gcc/13.3.0 python/3.11.9


TRAIT=$1
LDCTS=$2
OUT_NUM=$3
ANC=${4:-EUR}

case "$ANC" in
  EUR|EAS) ;;
  *) echo "Ancestry must be EUR or EAS" >&2; exit 1 ;;
esac

SCRIPT_PATH=$(dirname $0)
DATA_PATH=$SCRIPT_PATH/../data

echo $DATA_PATH

#SLDSC="python /home1/artemkim/CTFM_3_test/install/ldsc-2.0.1/ldsc.py";
SLDSC="ldsc.py";
REF_DIR="$DATA_PATH/1000G_${ANC}_Phase3"
BASELINE="$REF_DIR/baseline_v1.2/baseline.";
WEIGHTS="--w-ld-chr $REF_DIR/weights/weights.hm3_noMHC.";


# modify SUMSTAT and OUTDIR variables as needed
SUMSTAT="sumstats_ready/$TRAIT.sumstats.gz"
OUTDIR="out/SLDSC/$TRAIT/"
if [ "$ANC" = "EAS" ]; then OUTDIR="out/SLDSC_EAS/$TRAIT/"; fi

OUTFILE="$OUTDIR/$TRAIT.$OUT_NUM"
mkdir -p $OUTDIR

#sed "s@data@${DATA_PATH}@g" $SCRIPT_PATH/../$LDCTS > $OUTFILE.cts

sed "s@data@${DATA_PATH}@g" "$LDCTS" > "$OUTFILE.cts" || exit $?


$SLDSC \
  --h2-cts $SUMSTAT \
  --ref-ld-chr $BASELINE \
  --ref-ld-chr-cts $OUTFILE.cts \
  $WEIGHTS \
  --out $OUTFILE || exit $?

rm $OUTFILE.cts

