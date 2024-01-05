ARG BASE_IMAGE
FROM ${BASE_IMAGE}

ENV GP_EXT_DEPS=/gridpack-dependencies
ENV BOOST_VERSION "1.78.0"
ENV GA_VERSION "5.8"
ENV PETSC_VERSION "3.16.4"
ENV LD_LIBRARY_PATH=${GP_EXT_DEPS}/boost/install_for_gridpack/lib:${GP_EXT_DEPS}/ga/install_for_gridpack/lib:${GP_EXT_DEPS}/petsc/install_for_gridpack/lib:${LD_LIBRARY_PATH}

WORKDIR ${GP_EXT_DEPS}

COPY *.sh .

RUN (./install_environment_packages.sh && ./install_gridpack_deps.sh) || find ${GP_EXT_DEPS} -type f -regex ".*log.*" && rm *.sh
