#!/bin/bash

echo "Build and install annotateCSF"
echo "============================="

echo "Run with '$0 -create' to build installer package"

if [ ! -f "files/ref2.h5ad" ]; then
    echo "Preparing data"
    echo "--------------"
    wget 'https://zenodo.org/record/6912278/files/matrix.mtx?download=1' -O 'test_data/10x_csf/matrix.mtx' 
    gzip -d files/ref2.h5ad.gz
fi

if [ ! -d "env" ]; then
    echo "Installing required packages"
    echo "----------------------------"
    sudo add-apt-repository ppa:deadsnakes/ppa
    sudo apt-get update
    sudo apt-get install python3.9 python3.9-tk python3.9-venv python3.9-dev pip  # virtualenv

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
    pyinstaller --clean -y -n "run_aCSF" --copy-metadata anndata --copy-metadata scvi-tools --copy-metadata scanpy --hidden-import PIL --hidden-import PIL._imagingtk --hidden-import PIL._tkinter_finder --add-data="env/lib/python3.9/site-packages/numpy.libs/*:numpy.libs" --add-data="env/lib/python3.9/site-packages/scipy.libs/*:scipy.libs" --add-data="env/lib/python3.9/site-packages/Pillow.libs/*:Pillow.libs" --add-data="env/lib/python3.9/site-packages/h5py.libs/*:h5py.libs" --add-data="files/*:files" annotateCSF.py

    echo "Launching distributable version"
    echo "-------------------------------"
    echo "Waiting for annotateCSF to be closed . . . "
    ./dist/run_aCSF/run_aCSF
fi

echo "Finished."
