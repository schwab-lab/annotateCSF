echo "Build and install annotateCSF"
echo "============================="

echo "Installing required packages"
echo "----------------------------"
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt-get update
sudo apt-get install python3.9 python3.9-tk python3.9-venv python3.9-dev pip  # virtualenv

echo "Preparing data"
echo "--------------"
gzip -d files/ref2.h5ad.gz
wget 'https://zenodo.org/record/6795112/files/matrix.mtx?download=1' -O 'test_data\10x_csf\matrix.mtx' 

echo "Creating virtul environment"
echo "---------------------------"
#virtualenv --python="/usr/bin/python3.9" env
python3.9 -m venv env
source env/bin/activate

echo "Installing dependencies"
echo "-----------------------"
pip install --no-cache-dir -r requirements_linux.txt
pip install --no-cache-dir torch==1.9.1+cu111

echo "Launching annotateCSF"
echo "---------------------"
echo "The tool will be launched now. Please check if it works before the script continues!"
echo "Waiting for annotateCSF to be closed . . . "
python annotateCSF.py

echo "Installing pyInstaller"
echo "----------------------"
pip --no-cache-dir install pyinstaller==5.0.1

echo "Compiling distributable version"
echo "-------------------------------"
pyinstaller --clean -y -n "run_aCSF" --copy-metadata anndata --copy-metadata scvi-tools --copy-metadata scanpy annotateCSF.py

echo "Launching distributable version"
echo "-------------------------------"
echo "Waiting for annotateCSF to be closed . . . "
./dist/run_aCSF/run_aCSF

echo "Finished."
