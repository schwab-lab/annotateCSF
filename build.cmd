@echo off
cd /d %~dp0
python -m venv env
env\Scripts\activate.bat
pip install -r requirements_win.txt
pip install torch==1.9.1+cu111
python annotateCSF.py

pause

pip install pyinstaller==5.0.1
rem pyinstaller --clean -y -n "run_aCSF" --hidden-import="sklearn.utils._typedefs" --hidden-import="sklearn.utils._heap" --hidden-import="sklearn.utils._sorting" --hidden-import="sklearn.metrics._pairwise_distances_reduction" --hidden-import="sklearn.utils._vector_sentinel" --hidden-import="sklearn.utils._cython_blas" --hidden-import="sklearn.neighbors.typedefs" --hidden-import="sklearn.neighbors.quad_tree" --hidden-import="sklearn.tree._utils" --hidden-import="sklearn.neighbors._typedefs"  --hidden-import="sklearn.neighbors._partition_nodes" --add-data="files/acsf_logo.png;files" --copy-metadata anndata --copy-metadata scvi-tools  --copy-metadata scanpy --add-data="files/ref2.h5ad;files" --add-data="files/gene_symbols_to_match.csv;files" --icon="files/icon3.ico" --add-data="files/icon3.ico;files" --add-data="files/plot_options.png;files"  --add-data="files/dge_options.png;files" --add-data="files/acsf_logo_basic_workflow.png;files" annotateCSF.py
pyinstaller --clean -F -y -n "run_aCSF" --copy-metadata anndata --copy-metadata scvi-tools --copy-metadata scanpy --add-data="env/Lib/site-packages/tables.libs/*;./tables.libs" --add-data="files/*;files" --icon="files/icon3.ico" annotateCSF.py

env\Scripts\deactivate.bat

dist\run_aCSF\run_aCSF.exe

pause