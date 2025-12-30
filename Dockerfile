FROM ubuntu:questing

ENV DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC

RUN apt-get update && \
    apt-get install -y --no-install-recommends make wget tzdata git \
    python3 python3-pip python3-venv \
    build-essential openmpi-bin openmpi-common openmpi-doc libopenmpi-dev && \
    apt-get clean

COPY README.md install_gridpack.sh install_gridpack_deps.sh /app/
COPY docs /app/docs
COPY python /app/python
COPY src /app/src
