SCRIPT_PATH=$(dirname $0)
DATA_PATH=$SCRIPT_PATH/../data
cd $DATA_PATH

mkdir -p SLDSC_default/
cd SLDSC_default

echo "Downloading S-LDSC bed files CT-FM-SNP"


wget https://zenodo.org/records/11200933/files/SLDSC_beds.zip?download=1  -O SLDSC_beds.zip
unzip SLDSC_beds.zip

rm -r __MACOSX # need to update archives and remove MACOSX junk files


echo "Done"