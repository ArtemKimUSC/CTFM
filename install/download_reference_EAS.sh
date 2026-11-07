
SCRIPT_PATH=$(dirname $0)
DATA_PATH=$SCRIPT_PATH/../data



cd $DATA_PATH

echo "Downloading EAS 1000G Reference files for CT-FM"

wget https://zenodo.org/records/23923165/files/1000G_EAS_Phase3.tar.gz?download=1 -O CTFM_EAS_1000G.tar.gz

tar -zxvf CTFM_EAS_1000G.tar.gz

rm CTFM_EAS_1000G.tar.gz

echo "Done"



echo "Downloading EAS S-LDSC reference annotations for CT-FM"

wget https://zenodo.org/records/24018537/files/SLDSC_EAS.tar.gz?download=1 -O SLDSC_EAS.tar.gz
tar -zxvf SLDSC_EAS.tar.gz

rm SLDSC_EAS.tar.gz

echo "Done"

rm -r __MACOSX # need to update archives and remove MACOSX junk files
