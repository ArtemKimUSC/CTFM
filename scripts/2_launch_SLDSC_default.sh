TRAIT=$1
ANC=${2:-EUR}

case "$ANC" in
  EUR) PANEL=SLDSC_default ;;
  EAS) PANEL=SLDSC_EAS ;;
  *) echo "Ancestry must be EUR or EAS" >&2; exit 1 ;;
esac

SCRIPT_PATH=$(dirname $0)
DATA_PATH=$SCRIPT_PATH/../data

bash "$SCRIPT_PATH/sldsc.sh" "$TRAIT" "$DATA_PATH/$PANEL/sldsc_annots_00" 0 "$ANC" || exit $?
bash "$SCRIPT_PATH/sldsc.sh" "$TRAIT" "$DATA_PATH/$PANEL/sldsc_annots_01" 1 "$ANC" || exit $?
bash "$SCRIPT_PATH/sldsc.sh" "$TRAIT" "$DATA_PATH/$PANEL/sldsc_annots_02" 2 "$ANC" || exit $?
bash "$SCRIPT_PATH/sldsc.sh" "$TRAIT" "$DATA_PATH/$PANEL/sldsc_annots_03" 3 "$ANC" || exit $?
bash "$SCRIPT_PATH/sldsc.sh" "$TRAIT" "$DATA_PATH/$PANEL/sldsc_annots_04" 4 "$ANC" || exit $?

