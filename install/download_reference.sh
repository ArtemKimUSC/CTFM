
SCRIPT_PATH=$(dirname $0)
DATA_PATH=$SCRIPT_PATH/../data



cd $DATA_PATH

echo "Downloading 1000G Reference files for CT-FM"

wget https://zenodo.org/records/11193775/files/CTFM_1000G.zip?download=1 -O CTFM_1000G.zip

unzip CTFM_1000G.zip

rm CTFM_1000G.zip

echo "Done"



echo "Downloading S-LDSC reference annotations for CT-FM"

wget https://zenodo.org/records/12192853/files/SLDSC_default.zip?download=1  -O SLDSC_default.zip
unzip SLDSC_default.zip

rm SLDSC_default.zip

echo "Done"

rm -r __MACOSX # need to update archives and remove MACOSX junk files
