ARG BASE_IMAGE
FROM ${BASE_IMAGE}

ENV GP_EXT_DEPS=/gridpack-dependencies
ENV BOOST_VERSION="1.78.0"
ENV GA_VERSION="5.8"
ENV PETSC_VERSION="3.16.4"
ENV LD_LIBRARY_PATH=${GP_EXT_DEPS}/boost/install_for_gridpack/lib:${GP_EXT_DEPS}/ga/install_for_gridpack/lib:${GP_EXT_DEPS}/petsc/install_for_gridpack/lib:${LD_LIBRARY_PATH}

ARG http_proxy
ARG https_proxy
ENV http_proxy=${http_proxy}
ENV https_proxy=${https_proxy}

WORKDIR ${GP_EXT_DEPS}

COPY install_environment_packages.sh .
RUN ./install_environment_packages.sh && rm *.sh

COPY install_gridpack_deps.sh .
RUN ./install_gridpack_deps.sh && rm *.sh
