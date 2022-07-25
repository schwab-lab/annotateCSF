In case, you want to manually install annotateCSF in your system (NOT recommended!), this are the minimal requirements.
Please refer to Readme.md for recommended installation options.


# Python & library release versions (Windows)

#python
python3.9
python3.9-tk
python3.9-dev
python3.9-venv
python-is-python3

#libraries
pandas==1.3.3
pandastable==0.13.0
numpy==1.21.5
scipy==1.7.1
statsmodels==0.13.0rc0
matplotlib==3.4.3
seaborn==0.11.2
scrublet==0.2.3
PIL==8.3.2
torch (& cuda) == 1.9.1+cu111
scvi (scvi-tools) == 0.13.0
scanpy==1.8.1
skmisc==0.1.4
termcolor==1.1.0

#pyinstaller windows
pyinstaller==5.0.1
 
# code pyinstaller:
pyinstaller --clean -F -y -n "run_aCSF" --copy-metadata anndata --copy-metadata scvi-tools --copy-metadata scanpy --add-data="env/Lib/site-packages/tables.libs/*;./tables.libs" --add-data="files/*;files" --icon="files/icon3.ico" annotateCSF.py