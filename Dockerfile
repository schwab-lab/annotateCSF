FROM debian:bullseye-slim

## Arguments for sharing X-server
## (from: https://riptutorial.com/docker/example/21831/running-gui-apps-in-a-linux-container)
#ARG user
#ARG uid
#ARG gid

RUN apt update && apt install -y python3.9 python3.9-tk python3.9-dev gcc g++ libblas-dev liblapack-dev wget

WORKDIR /annotatecsf

COPY requirements_linux.txt requirements_linux.txt

RUN pip install --no-cache-dir -r requirements_linux.txt
#RUN pip install --no-cache-dir torch==1.9.1+cu111

COPY files/* files/
COPY annotateCSF.py annotateCSF.py


RUN echo '#!/bin/bash' > build.sh \
 && echo 'pip --no-cache-dir install pyinstaller==5.0.1' >> build.sh \
 && echo '#wget https://raw.githubusercontent.com/uni-ms/annotateCSF/main/annotateCSF.py' >> build.sh \
 && echo 'pyinstaller --clean -y -n "run_aCSF" --copy-metadata anndata --copy-metadata scvi-tools --copy-metadata scanpy --hidden-import PIL --hidden-import PIL._imagingtk --hidden-import PIL._tkinter_finder --add-data=".local/lib/python3.9/site-packages/numpy.libs/*:." --add-data=".local/lib/python3.9/site-packages/scipy.libs/*:." --add-data=".local/lib/python3.9/site-packages/Pillow.libs/*:." --add-data=".local/lib/python3.9/site-packages/h5py.libs/*:." --add-data="files/*:." annotateCSF.py' >> build.sh

RUN echo '#!/bin/bash' > run_user.sh \
 && echo 'useradd -m $USERNAME' >> run_user.sh \
 && echo '$USERNAME:$USERNAME | chpasswd' >> run_user.sh \
 && echo 'usermod --shell /bin/bash $USERNAME' >> run_user.sh \
 && echo 'usermod --uid ${uid} $USERNAME' >> run_user.sh \
 && echo 'groupmod --gid ${gid} $USERNAME' >> run_user.sh \
 && echo 'chown $USERNAME:$USERNAME /annotatecsf' >> run_user.sh \
 && echo 'su $USERNAME' >> run_user.sh \
 && echo 'cd /annotatecsf' >> run_user.sh \
 && echo '#wget https://raw.githubusercontent.com/uni-ms/annotateCSF/main/annotateCSF.py' >> run_user.sh \
 && echo 'python annotateCSF.py' >> run_user.sh


CMD python annotateCSF.py
