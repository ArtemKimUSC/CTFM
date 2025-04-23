
SCRIPT_PATH=$(dirname $0)
DATA_PATH=$SCRIPT_PATH/../data



cd $DATA_PATH

echo "Downloading S-LDSC cancer annotations for CT-FM"

wget https://zenodo.org/records/15271630/files/SLDSC_cancer.tar.gz?download=1 -O SLDSC_cancer.tar.gz

tar -xvzf SLDSC_cancer.tar.gz

rm SLDSC_cancer.tar.gz

echo "Done"
