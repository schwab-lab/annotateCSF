# WINDOWS Python & library release versions

#python
Python3.9.7

#pyinstaller windows
pyinstaller == 5.0.1

#libraries
tkinter == 8.6
PIL == 8.3.2
pandas == 1.3.3
scvi (scvi-tools) == 0.13.0
scanpy == 1.8.1
skmisc == 0.1.4
matplotlib == 3.4.3
seaborn == 0.11.2
scrublet == 0.2.3
torch (& cuda) == 1.9.1+cu111
numpy == 1.21.5
statsmodels == 0.13.0rc0
pandastable == 0.13.0
termcolor == 1.1.0
scipy == 1.7.1
 
# code pyinstaller:
!! Pfade müssen angepasst werden
pyinstaller --clean -y -n "run_aCSF" --add-data="E:/py_guis/tkinter/acsf_logo.png;files" --copy-metadata anndata --copy-metadata scvi-tools  --copy-metadata scanpy --add-data="E:/py_guis/tkinter/adata_ref/ref2.h5ad;files" --add-data="E:/py_guis/tkinter/gene_symbols_to_match.csv;files" --icon="E:/py_guis/tkinter/icon3.ico" --add-data="E:/py_guis/tkinter/icon3.ico;files" --add-data="E:/py_guis/tkinter/plot_options.png;files"  --add-data="E:/py_guis/tkinter/dge_options.png;files" --add-data="E:/py_guis/tkinter/acsf_logo_basic_workflow.png;files" E:/py_guis/tkinter/annotateCSF.py

