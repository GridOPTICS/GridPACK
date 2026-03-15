Docker
======

The GridPACK CLI operates inside Docker containers without modification.
This chapter describes how to build a container image with CLI support
and how to run simulations within a containerized environment.

Running the CLI in Docker
-------------------------

Once GridPACK is built and installed inside a Docker image, the
``gridpack`` command is available on the container's PATH.

Basic invocation::

    docker run gridpack-image gridpack powerflow /data/input.xml

Mount a host directory to provide input files and collect output::

    docker run -v $(pwd)/data:/data gridpack-image \
        gridpack pf /data/input.xml --output /data/results.txt

Parallel execution with MPI inside the container::

    docker run gridpack-image \
        mpiexec -np 4 gridpack dsf /data/input.xml

Building a Docker Image
------------------------

The following Dockerfile demonstrates how to build a GridPACK image
with Python CLI support:

.. code-block:: dockerfile

    FROM ubuntu:22.04

    # Install system dependencies
    RUN apt-get update && apt-get install -y \
        build-essential cmake \
        libopenmpi-dev openmpi-bin \
        python3 python3-pip python3-dev \
        libboost-all-dev libpetsc-dev \
        && rm -rf /var/lib/apt/lists/*

    RUN pip3 install mpi4py

    # Build GridPACK with Python bindings
    COPY . /gridpack
    WORKDIR /gridpack/src/build
    RUN cmake .. -DENABLE_PYTHON=ON \
        && make -j$(nproc) \
        && make install

    # Install the Python CLI package
    WORKDIR /gridpack/python
    RUN pip3 install .

    # Verify
    RUN gridpack --version

Build and run::

    docker build -t gridpack .
    docker run -v $(pwd)/data:/data gridpack powerflow /data/input.xml

Practical Considerations
------------------------

- **Input files.** XML configuration files, PSS/E RAW files, DYR files,
  and measurement files must be accessible inside the container. Use
  volume mounts (``-v``) to map host directories into the container
  filesystem.

- **Output files.** Files written by the application (power flow
  results, contingency output files, generator watch data) are created
  inside the container. Mount an output directory to persist results
  on the host.

- **MPI parallelism.** The number of MPI processes is controlled by
  ``mpiexec`` inside the container, not by Docker scaling or
  replicas. For multi-node MPI execution, use Docker Compose or a
  container orchestration platform with appropriate MPI networking.

- **Python version.** The compiled Python extension (``.so``) must
  match the Python interpreter version inside the container.
  Ensure consistency between the build-time and runtime Python
  installations.
