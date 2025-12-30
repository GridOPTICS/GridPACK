FROM ubuntu:questing

ARG boost_version=1.81.0
ARG ga_version=5.9.1
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

