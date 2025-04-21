TRAIT=$1

SCRIPT_PATH=$(dirname $0)
DATA_PATH=$SCRIPT_PATH/../data

sbatch --mem=40G -t 0-24:00:00 --wrap="bash scripts/sldsc.sh $TRAIT $DATA_PATH/SLDSC_default/sldsc_annots_00 0"
sbatch --mem=40G -t 0-24:00:00 --wrap="bash scripts/sldsc.sh $TRAIT $DATA_PATH/SLDSC_default/sldsc_annots_01 1"
sbatch --mem=40G -t 0-24:00:00 --wrap="bash scripts/sldsc.sh $TRAIT $DATA_PATH/SLDSC_default/sldsc_annots_02 2"
sbatch --mem=40G -t 0-24:00:00 --wrap="bash scripts/sldsc.sh $TRAIT $DATA_PATH/SLDSC_default/sldsc_annots_03 3"
sbatch --mem=40G -t 0-24:00:00 --wrap="bash scripts/sldsc.sh $TRAIT $DATA_PATH/SLDSC_default/sldsc_annots_04 4"
