FROM python:3.9-slim-bullseye

# Arguments for sharing X-server
# (from: https://riptutorial.com/docker/example/21831/running-gui-apps-in-a-linux-container)
ARG user
ARG uid
ARG gid

RUN apt update && apt install -y python3.9-tk python3.9-dev gcc g++ libblas-dev liblapack-dev

#Add new user with our credentials
ENV USERNAME ${user}
RUN useradd -m $USERNAME && \
      echo "$USERNAME:$USERNAME" | chpasswd && \
      usermod --shell /bin/bash $USERNAME && \
      usermod  --uid ${uid} $USERNAME && \
      groupmod --gid ${gid} $USERNAME

USER ${user}

WORKDIR /home/${user}


COPY annotateCSF.py annotateCSF.py
COPY requirements_linux.txt requirements_linux.txt
COPY files/* files/

RUN pip install --no-cache-dir -r requirements_linux.txt
#RUN pip install --no-cache-dir torch==1.9.1+cu111

CMD python annotateCSF.py

