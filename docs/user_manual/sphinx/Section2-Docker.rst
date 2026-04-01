Using GridPACK with Docker
===========================

GridPACK is available as a Docker image with all dependencies pre-installed and ready to use.

Quick Start
-----------

Pull the Image
~~~~~~~~~~~~~~

.. code-block:: bash

   docker pull pnnl/gridpack:latest

Run with Your Files
~~~~~~~~~~~~~~~~~~~

The container starts in ``/app/workspace`` by default. Mount your current directory to work with your files:

.. code-block:: bash

   docker run -it --rm -v $(pwd):/app/workspace pnnl/gridpack:latest

Usage Examples
--------------

Interactive Shell
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Mount your working directory and start an interactive session
   docker run -it --rm -v $(pwd):/app/workspace pnnl/gridpack:latest bash

Run a Specific Command
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Run a GridPACK application directly
   docker run --rm -v $(pwd):/app/workspace pnnl/gridpack:latest \
     pf.x input_14.xml

Python Scripts
~~~~~~~~~~~~~~

.. code-block:: bash

   # Execute a Python script using GridPACK Python bindings
   docker run --rm -v $(pwd):/app/workspace pnnl/gridpack:latest \
     python3 my_gridpack_script.py

Run Tests
~~~~~~~~~

.. code-block:: bash

   # Run the GridPACK test suite
   docker run --rm pnnl/gridpack:latest bash -c "cd /app/src/build && ctest"

Working Directory
-----------------

The container is configured with:

- **Working directory**: ``/app/workspace``
- **GridPACK installation**: ``/app/src/install``
- **Build directory**: ``/app/src/build``

Mount your data/scripts to ``/app/workspace`` to access them from within the container.

Environment Variables
---------------------

The following environment variables are pre-configured:

.. code-block:: bash

   GRIDPACK_DIR=/app/src/install
   PATH=/app/src/install/bin:$PATH
   LD_LIBRARY_PATH=/deps/boost-1.81.0/install_for_gridpack/lib:/deps/ga-5.9.1/install_for_gridpack/lib:/deps/petsc/install_for_gridpack/lib

Multi-Architecture Support
---------------------------

The Docker image supports both AMD64 and ARM64 architectures. Docker automatically pulls the correct image for your platform.

To explicitly specify an architecture:

.. code-block:: bash

   # For AMD64 (Intel/AMD processors)
   docker pull --platform linux/amd64 pnnl/gridpack:latest

   # For ARM64 (Apple Silicon, AWS Graviton)
   docker pull --platform linux/arm64 pnnl/gridpack:latest

Available Tags
--------------

- ``latest`` - Latest stable release
- ``dev`` - Development builds
- ``vX.Y.Z`` - Specific version releases

Tips
----

- Add ``-it`` for interactive terminal access
- Add ``--rm`` to automatically remove the container after exit
- Use ``-v`` to mount local directories into the container
- Increase memory limits for large simulations: ``--memory=8g``
- Use multiple cores: Set ``OMPI_NUM_PROCS`` or pass to MPI commands

Example Workflow
----------------

.. code-block:: bash

   # 1. Create a project directory
   mkdir my-gridpack-project
   cd my-gridpack-project

   # 2. Create your input files
   cat > input.xml << EOF
   <?xml version="1.0" encoding="utf-8"?>
   <Configuration>
     <!-- Your GridPACK configuration -->
   </Configuration>
   EOF

   # 3. Run GridPACK in the container
   docker run -it --rm -v $(pwd):/app/workspace pnnl/gridpack:latest bash

   # Inside container:
   # powerflow.x input.xml
   # exit

Troubleshooting
---------------

Permission Issues
~~~~~~~~~~~~~~~~~

If you encounter permission issues with mounted volumes:

.. code-block:: bash

   # Run with your user ID
   docker run -it --rm -u $(id -u):$(id -g) -v $(pwd):/app/workspace pnnl/gridpack:latest

MPI Warnings
~~~~~~~~~~~~

The container sets ``OMPI_ALLOW_RUN_AS_ROOT=1`` to allow running MPI as root. This is safe for containerized environments.

Running Applications with MPI
------------------------------

GridPACK applications can be run with MPI inside the Docker container. For example:

.. code-block:: bash

   # Run with 4 MPI processes
   docker run --rm -v $(pwd):/app/workspace pnnl/gridpack:latest \
     mpirun -n 4 powerflow.x input.xml

For detailed information on MPI execution, application arguments, and input file conventions,
see :doc:`Section3-ConfigureBuild` (Running GridPACK Applications section).

Building from Source
--------------------

If you need to build GridPACK from source instead of using Docker, see :doc:`Section3-ConfigureBuild`
for detailed instructions on using build scripts or manual compilation.
