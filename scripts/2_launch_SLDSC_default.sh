TRAIT=$1

SCRIPT_PATH=$(dirname $0)
DATA_PATH=$SCRIPT_PATH/../data

bash scripts/sldsc.sh $TRAIT $DATA_PATH/SLDSC_default/sldsc_annots_00 0
bash scripts/sldsc.sh $TRAIT $DATA_PATH/SLDSC_default/sldsc_annots_01 1
bash scripts/sldsc.sh $TRAIT $DATA_PATH/SLDSC_default/sldsc_annots_02 2
bash scripts/sldsc.sh $TRAIT $DATA_PATH/SLDSC_default/sldsc_annots_03 3
bash scripts/sldsc.sh $TRAIT $DATA_PATH/SLDSC_default/sldsc_annots_04 4
