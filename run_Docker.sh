# Build docker container (optional)
docker build --build-arg user=$USER --build-arg uid=$(id -u) --build-arg gid=$(id -g) -t unims/annotatecsf:1.0 -f Dockerfile .

# Use Docker container to create the installer
mkdir output
docker run -v ${PWD}/output:/home/$USER/dist unims/annotatecsf:1.0 bash -c 'pip --no-cache-dir install pyinstaller==5.0.1 && ~/.local/bin/pyinstaller --clean -y -n "run_aCSF" --copy-metadata anndata --copy-metadata scvi-tools --copy-metadata scanpy --hidden-import PIL --hidden-import PIL._imagingtk --hidden-import PIL._tkinter_finder --add-data=".local/lib/python3.9/site-packages/numpy.libs/*:numpy.libs" --add-data=".local/lib/python3.9/site-packages/scipy.libs/*:scipy.libs" --add-data=".local/lib/python3.9/site-packages/Pillow.libs/*:Pillow.libs" --add-data=".local/lib/python3.9/site-packages/h5py.libs/*:h5py.libs" --add-data="files/*:files" annotateCSF.py'

# Simple, but insecure method
#xhost +local:root

# Now, before spawning a GUI container, we have to create a xauth file with access permission
xauth nlist $DISPLAY | sed -e 's/^..../ffff/' | xauth -f /tmp/.docker.xauth nmerge -

# For Windows: DISPLAY=host.docker.internal:0
# For MacOS: DISPLAY=docker.for.mac.host.internal:0

# This file has to be mounted into the container when creating/running it
docker run -it -v ${PWD}/test_data:/home/$USER/test_data -e DISPLAY=unix$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v /tmp/.docker.xauth:/tmp/.docker.xauth:rw -e XAUTHORITY=/tmp/.docker.xauth unims/annotatecsf:1.0

#xhost -local:root

