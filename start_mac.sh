#!/bin/bash

echo "Build and install annotateCSF"
echo "============================="

if [ ! -f "files/ref2.h5ad" ]; then
    echo "Preparing data"
    echo "--------------"
    wget 'https://zenodo.org/record/6912278/files/matrix.mtx?download=1' -O 'test_data/10x_csf/matrix.mtx' 
    gzip -d files/ref2.h5ad.gz
fi

if [ ! -d "env" ]; then
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
    pip install torch==1.9.1+cu111 -f https://download.pytorch.org/whl/torch_stable.html
fi

source env/bin/activate

echo "Launching annotateCSF"
echo "---------------------"
echo "The tool will be launched now. Please check if it works before the script continues!"
echo "Waiting for annotateCSF to be closed . . . "
python annotateCSF.py


if [ "$1" == "-create" ]; then
    echo "Installing pyInstaller"
    echo "----------------------"
    pip install --no-cache-dir pyinstaller==5.0.1

    echo "Compiling distributable version"
    echo "-------------------------------"
    pyinstaller --clean -y -n "run_aCSF" --target-arch x86_64 --copy-metadata anndata --copy-metadata scvi-tools --copy-metadata scanpy --add-data="files/*:files" --icon="files/icon3.ico" annotateCSF.py

    echo "Launching distributable version"
    echo "-------------------------------"
    echo "Waiting for annotateCSF to be closed . . . "
    ./dist/run_aCSF/run_aCSF
fi

echo "Finished."
