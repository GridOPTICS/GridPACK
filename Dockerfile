FROM ubuntu:questing

ARG boost_version=1.81.0
ARG ga_version=5.9.1
ARG petsc_version=3.24.2

ENV DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC GNUMAKEFLAGS=--no-print-directory
ENV PETSC_DIR=/deps/petsc PETSC_ARCH=build-dir
ENV OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# GridPACK dependency environment variables
ENV GRIDPACK_ROOT_DIR=/app
ENV GP_EXT_DEPS=/deps
ENV LD_LIBRARY_PATH=/deps/boost-${boost_version}/install_for_gridpack/lib:/deps/ga-${ga_version}/install_for_gridpack/lib:/deps/petsc/install_for_gridpack/lib
ENV DYLD_LIBRARY_PATH=/deps/boost-${boost_version}/install_for_gridpack/lib:/deps/ga-${ga_version}/install_for_gridpack/lib:/deps/petsc/install_for_gridpack/lib

RUN apt-get update && \
    apt-get install -y --no-install-recommends cmake make wget tzdata git gfortran \
    python3 python3-pip python3-venv \
    build-essential openmpi-bin openmpi-common openmpi-doc libopenmpi-dev && \
    apt-get clean

COPY README.md install_gridpack.sh /app/
COPY docs /app/docs
COPY python /app/python
COPY src /app/src

# Compile/Install Boost
WORKDIR /deps
RUN wget "https://github.com/boostorg/boost/releases/download/boost-${boost_version}/boost-${boost_version}.tar.gz"
RUN tar -xf "boost-${boost_version}.tar.gz"
WORKDIR "/deps/boost-${boost_version}"
RUN ./bootstrap.sh --prefix=install_for_gridpack --with-libraries=mpi,serialization,random,filesystem,system
RUN echo 'using mpi : mpicxx ; ' >> project-config.jam
RUN ./b2 -a -d+2 link="shared" stage
RUN ./b2 -a -d+2 link="shared" install

# Compile/Install GA (Global Arrays)
WORKDIR /deps
RUN wget "https://github.com/GlobalArrays/ga/releases/download/v${ga_version}/ga-${ga_version}.tar.gz"
RUN tar -xf "ga-${ga_version}.tar.gz"
WORKDIR "/deps/ga-${ga_version}"
RUN ./configure --with-mpi-ts --disable-f77 \
    --without-blas --without-lapack --without-scalapack \
    --enable-cxx --enable-i4 \
    --prefix=${PWD}/install_for_gridpack \
    CFLAGS="-Wno-implicit-function-declaration -Wno-incompatible-pointer-types -Wno-old-style-definition" \
    CXXFLAGS="-Wno-incompatible-pointer-types" \
    --enable-shared=yes --enable-static=no
RUN make -j 10 install
WORKDIR /deps

# Compile/Install PETSc
RUN git clone https://gitlab.com/petsc/petsc.git
WORKDIR /deps/petsc
RUN git checkout "tags/v${petsc_version}" -b "v${petsc_version}"
RUN ./configure \
    --prefix=${PWD}/install_for_gridpack \
    --scalar-type=real  \
    --with-fortran-bindings=0 \
    --download-superlu_dist \
    --download-metis \
    --download-parmetis \
    --download-suitesparse \
    --download-f2cblaslapack \
    --download-scalapack \
    --download-mumps \
    --download-cmake=0 \
    --with-sowing=0 \
    --with-debugging=0 \
    --with-shared-libraries=1
RUN make all
RUN make install
RUN make PETSC_DIR=/deps/petsc/install_for_gridpack PETSC_ARCH="" check

WORKDIR /app