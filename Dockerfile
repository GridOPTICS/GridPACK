FROM ubuntu:questing

ARG boost_version=1.81.0
ARG ga_version=5.9.1
ARG petsc_version=3.24.2

ENV DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC GNUMAKEFLAGS=--no-print-directory
ENV OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

ENV GRIDPACK_ROOT_DIR=/app GP_EXT_DEPS=/deps
ENV GRIDPACK_INSTALL_DIR=${GRIDPACK_ROOT_DIR}/src/install GRIDPACK_BUILD_DIR=${GRIDPACK_ROOT_DIR}/src/build
ENV GRIDPACK_DIR=${GRIDPACK_INSTALL_DIR}

ENV boost_dir=${GP_EXT_DEPS}/boost-${boost_version} \
    ga_dir=${GP_EXT_DEPS}/ga-${ga_version} \
    petsc_dir=${GP_EXT_DEPS}/petsc

ENV boost_gp_dir=${boost_dir}/install_for_gridpack \
    ga_gp_dir=${ga_dir}/install_for_gridpack \
    petsc_gp_dir=${petsc_dir}/install_for_gridpack

ENV PETSC_DIR=${petsc_dir} PETSC_ARCH=build-dir

ENV LD_LIBRARY_PATH=${boost_gp_dir}/lib:${ga_gp_dir}/lib:${petsc_gp_dir}/lib
ENV DYLD_LIBRARY_PATH=${LD_LIBRARY_PATH}

RUN apt-get update && \
    apt-get install -y --no-install-recommends cmake make wget tzdata git gfortran \
    build-essential openmpi-bin openmpi-common openmpi-doc libopenmpi-dev && \
    python3 python3-pip python3-venv python-is-python3 \
    apt-get clean

COPY README.md install_gridpack.sh ${GRIDPACK_ROOT_DIR}/
COPY docs ${GRIDPACK_ROOT_DIR}/docs
COPY python ${GRIDPACK_ROOT_DIR}/python
COPY src ${GRIDPACK_ROOT_DIR}/src

# Compile/Install Boost
WORKDIR ${GP_EXT_DEPS}
RUN wget "https://github.com/boostorg/boost/releases/download/boost-${boost_version}/boost-${boost_version}.tar.gz"
RUN tar -xf "boost-${boost_version}.tar.gz"
WORKDIR ${boost_dir}
RUN ./bootstrap.sh --prefix=install_for_gridpack --with-libraries=mpi,serialization,random,filesystem,system
RUN echo 'using mpi : mpicxx ; ' >> project-config.jam
RUN ./b2 -a -d+2 link="shared" stage
RUN ./b2 -a -d+2 link="shared" install

# Compile/Install GA (Global Arrays)
WORKDIR ${GP_EXT_DEPS}
RUN wget "https://github.com/GlobalArrays/ga/releases/download/v${ga_version}/ga-${ga_version}.tar.gz"
RUN tar -xf "ga-${ga_version}.tar.gz"
WORKDIR ${ga_dir}
RUN ./configure --with-mpi-ts --disable-f77 \
    --without-blas --without-lapack --without-scalapack \
    --enable-cxx --enable-i4 \
    --prefix=${ga_gp_dir} \
    CFLAGS="-Wno-implicit-function-declaration -Wno-incompatible-pointer-types -Wno-old-style-definition" \
    CXXFLAGS="-Wno-incompatible-pointer-types" \
    --enable-shared=yes --enable-static=no
RUN make -j 10 install
WORKDIR ${GP_EXT_DEPS}

# Compile/Install PETSc
RUN git clone https://gitlab.com/petsc/petsc.git
WORKDIR ${petsc_dir}
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
RUN make PETSC_DIR=${petsc_gp_dir} PETSC_ARCH="" check

WORKDIR /app