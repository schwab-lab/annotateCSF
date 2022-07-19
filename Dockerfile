FROM python:3.9.13-bullseye
#FROM debian:bullseye-slim

## Arguments for sharing X-server
## (from: https://riptutorial.com/docker/example/21831/running-gui-apps-in-a-linux-container)
ARG user=dockeruser
#ARG uid
#ARG gid
ENV USERNAME ${user}

RUN apt update && apt install -y python3.9 python3.9-tk python3.9-dev python3-pip gcc g++ libblas-dev liblapack-dev wget

RUN useradd --create-home --user-group $USERNAME \
 && echo "$USERNAME:$USERNAME" | chpasswd \
 && usermod --shell /bin/bash $USERNAME

USER ${user}
WORKDIR /home/${user}

COPY requirements_linux.txt requirements_linux.txt

RUN pip install --no-cache-dir -r requirements_linux.txt
#RUN pip install --no-cache-dir torch==1.9.1+cu111

COPY files/* files/
COPY annotateCSF.py annotateCSF.py


RUN echo '#!/bin/bash' > build.sh \
 && echo 'pip --no-cache-dir install pyinstaller==5.0.1' >> build.sh \
 && echo '#wget https://raw.githubusercontent.com/uni-ms/annotateCSF/main/annotateCSF.py' >> build.sh \
 && echo 'pyinstaller --clean -y -n "run_aCSF" --copy-metadata anndata --copy-metadata scvi-tools --copy-metadata scanpy --hidden-import PIL --hidden-import PIL._imagingtk --hidden-import PIL._tkinter_finder --add-data="/usr/local/lib/python3.9/dist-packages/numpy.libs/*:." --add-data="/usr/local/lib/python3.9/dist-packages/scipy.libs/*:." --add-data="/usr/local/lib/python3.9/dist-packages/Pillow.libs/*:." --add-data="/usr/local/lib/python3.9/dist-packages/h5py.libs/*:." --add-data="files/*:." annotateCSF.py' >> build.sh

RUN echo '#!/bin/bash' > run_user.sh \
 && echo 'usermod --uid ${UID} $USERNAME' >> run_user.sh \
 && echo 'groupmod --gid ${GID} $USERNAME' >> run_user.sh \
 && echo 'chown -R $USERNAME:$USERNAME /home/$USERNAME'
 && echo '#wget https://raw.githubusercontent.com/uni-ms/annotateCSF/main/annotateCSF.py' >> run_user.sh \
 && echo 'su $USERNAME -c "python annotateCSF.py"' >> run_user.sh


CMD python annotateCSF.py
