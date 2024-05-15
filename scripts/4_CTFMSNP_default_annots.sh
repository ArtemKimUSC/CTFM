SCRIPT_PATH=$(dirname $0)
DATA_PATH=$SCRIPT_PATH/../data

beds="$DATA_PATH/SLDSC_default/beds/bedlist.txt"

sed "s@data@${DATA_PATH}@g" $beds > beds.tmp

list_RS=$1
OUTDIR=$2



while read line; do

	SNPname=$(echo $line | awk '{print $4}')
	mkdir -p $OUTDIR

	FILE=$OUTDIR/$SNPname.out

	if test -f "$FILE"; then
    echo "$FILE already exists."
    continue
	fi


	echo "$line" > $OUTDIR/tmp.$SNPname.bed

	echo "Mapping annotations for SNP: $SNPname ..." 

	while read line; do
	
	bed=$(realpath "$line")
	echo $(basename $bed)
	bedtools intersect -a $OUTDIR/tmp.$SNPname.bed -b $bed -wa -wb >> $OUTDIR/$SNPname.out

	done < beds.tmp
	
	echo "Done for SNP: $SNPname"

	rm $OUTDIR/tmp.$SNPname.bed

done < $list_RS

rm beds.tmp

echo "DONE MAPPING SNPs, outdir:"$OUTDIR
