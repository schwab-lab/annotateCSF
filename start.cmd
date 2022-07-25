@echo off
cd /d %~dp0

echo Build and install annotateCSF
echo =============================
echo (Please check that Python 3.9 is installed!)
echo.
echo Run with '%0 -create' to build installer package
echo.

if not exist "files/ref2.h5ad" (
  echo Preparing data
  echo --------------
  powershell -command " (New-Object Net.WebClient).DownloadFile('https://zenodo.org/record/6795112/files/matrix.mtx?download=1', 'test_data\10x_csf\matrix.mtx') "
  python unzip.py -d files/ref2.h5ad.gz
)

if not exist env\ (
  echo Creating virtul environment
  echo ---------------------------
  python -m venv env
  call env\Scripts\activate.bat

  echo Installing dependencies
  echo -----------------------
  pip install --no-cache-dir -r requirements_win.txt
  pip install --no-cache-dir torch==1.9.1+cu111
)

call env\Scripts\activate.bat

echo Launching annotateCSF
echo ---------------------
echo The tool will be launched now. Please check if it works before the script continues!
echo Waiting for annotateCSF to be closed . . . 
start /w python annotateCSF.py


if "%1" == "-create" (
  echo Installing pyInstaller
  echo ----------------------
  pip --no-cache-dir install pyinstaller==5.0.1

  echo Compiling distributable version
  echo -------------------------------
  rem pyinstaller --clean -y -n "run_aCSF" --hidden-import="sklearn.utils._typedefs" --hidden-import="sklearn.utils._heap" --hidden-import="sklearn.utils._sorting" --hidden-import="sklearn.metrics._pairwise_distances_reduction" --hidden-import="sklearn.utils._vector_sentinel" --hidden-import="sklearn.utils._cython_blas" --hidden-import="sklearn.neighbors.typedefs" --hidden-import="sklearn.neighbors.quad_tree" --hidden-import="sklearn.tree._utils" --hidden-import="sklearn.neighbors._typedefs"  --hidden-import="sklearn.neighbors._partition_nodes" --add-data="files/acsf_logo.png;files" --copy-metadata anndata --copy-metadata scvi-tools  --copy-metadata scanpy --add-data="files/ref2.h5ad;files" --add-data="files/gene_symbols_to_match.csv;files" --icon="files/icon3.ico" --add-data="files/icon3.ico;files" --add-data="files/plot_options.png;files"  --add-data="files/dge_options.png;files" --add-data="files/acsf_logo_basic_workflow.png;files" annotateCSF.py
  pyinstaller --clean -y -n "run_aCSF" --win-no-prefer-redirects --copy-metadata anndata --copy-metadata scvi-tools --copy-metadata scanpy --add-data="env/Lib/site-packages/tables.libs/*;./tables.libs" --add-data="files/*;files" --icon="files/icon3.ico" annotateCSF.py

  echo Launching distributable version
  echo -------------------------------
  echo Waiting for annotateCSF to be closed . . . 
  start /w dist\run_aCSF\run_aCSF.exe
)

echo Cleanup environment
echo -------------------
call env\Scripts\deactivate.bat

echo Finished.
pause
