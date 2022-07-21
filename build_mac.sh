echo "Build and install annotateCSF"
echo "============================="

echo "Installing required packages"
echo "----------------------------"
brew install python@3.9 python-tk@3.9

echo "Creating virtul environment"
echo "---------------------------"
#virtualenv --python="/usr/bin/python3.9" env
python3.9 -m venv env
source env/bin/activate

echo "Installing dependencies"
echo "-----------------------"
pip install --no-cache-dir -r requirements_linux.txt
#pip install --no-cache-dir torch==1.9.1+cu111

echo "Launching annotateCSF"
echo "---------------------"
echo "The tool will be launched now. Please check if it works before the script continues!"
echo "Waiting for annotateCSF to be closed . . . "
python annotateCSF.py

echo "Installing pyInstaller"
echo "----------------------"
pip install --no-cache-dir pyinstaller==5.0.1

echo "Compiling distributable version"
echo "-------------------------------"
pyinstaller --clean -y -n "run_aCSF" --copy-metadata anndata --copy-metadata scvi-tools --copy-metadata scanpy --add-data="files/*:files" annotateCSF.py

echo "Launching distributable version"
echo "-------------------------------"
echo "Waiting for annotateCSF to be closed . . . "
./dist/run_aCSF/run_aCSF

echo "Finished."
