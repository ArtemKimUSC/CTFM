TRAIT=$1

SCRIPT_PATH=$(dirname $0)
DATA_PATH=$SCRIPT_PATH/../data

bash scripts/sldsc_cancer.sh $TRAIT $DATA_PATH/SLDSC_cancer/sldsc_annots_cancer_00 0
bash scripts/sldsc_cancer.sh $TRAIT $DATA_PATH/SLDSC_cancer/sldsc_annots_cancer_01 1
bash scripts/sldsc_cancer.sh $TRAIT $DATA_PATH/SLDSC_cancer/sldsc_annots_cancer_02 2
bash scripts/sldsc_cancer.sh $TRAIT $DATA_PATH/SLDSC_cancer/sldsc_annots_cancer_03 3
bash scripts/sldsc_cancer.sh $TRAIT $DATA_PATH/SLDSC_cancer/sldsc_annots_cancer_04 4
bash scripts/sldsc_cancer.sh $TRAIT $DATA_PATH/SLDSC_cancer/sldsc_annots_cancer_05 5
bash scripts/sldsc_cancer.sh $TRAIT $DATA_PATH/SLDSC_cancer/sldsc_annots_cancer_06 6