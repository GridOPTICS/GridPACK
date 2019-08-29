# -------------------------------------------------------------
# file: BuildBoostExternalProject.cmake
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created October 12, 2018 by William A. Perkins
# Last Change: 2018-10-19 13:08:22 d3g096
# -------------------------------------------------------------

# -------------------------------------------------------------
# BuildBoostExternalProject
# 
# This downloads and builds the specified Boost component
# libraries. The following variables are set to mimic
# find_package(Boost) variables:
#
#   Boost_VERSION
#   Boost_FOUND
#   Boost_INCLUDE_DIRS
#   Boost_LIBRARY_DIRS
#   Boost_LIBRARIES
#   Boost_<C>_FOUND
# -------------------------------------------------------------
function(BuildBoostExternalProject BOOST_COMPONENTS)

  # some requirements depending on what is being built 

  cmake_policy(SET CMP0057 NEW)
  
  if(mpi IN_LIST BOOST_COMPONENTS)
    find_package(MPI REQUIRED)
  endif()
  
  if(python IN_LIST BOOST_COMPONENTS)
    find_package(PythonInterp REQUIRED)
    find_package(PythonLibs REQUIRED)
    message(STATUS "PYTHONLIBS_VERSION_STRING=${PYTHONLIBS_VERSION_STRING}")
  endif()

  # This is the default version. It may get changed below.

  set(Boost_VERSION 1.65.1)
  set(BOOST_MD5 41d7542ce40e171f3f7982aff008ff0d)

  # set(Boost_VERSION 1.60.0)
  # set(BOOST_MD5 65a840e1a0b13a558ff19eeb2c4f0cbe)

  set(BOOST_OPTS "")
  set(BOOST_CONFIG_OPTS "")

  if(BUILD_SHARED_LIBS)
    list(APPEND BOOST_OPTS "link=shared")
  else()
    list(APPEND BOOST_OPTS "link=static")
  endif()
  
  # Figure out the Boost toolset from the C++ compiler that CMake found
  
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU") 
    set(BOOST_TOOLSET gcc)
    # Older GNU compilers need an older Boost version. I know 1.54 works.
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8")
      if (Boost_VERSION VERSION_GREATER "1.54.0")
        set(Boost_VERSION 1.54.0)
        set(BOOST_MD5 15cb8c0803064faef0c4ddf5bc5ca279)
        if(python IN_LIST BOOST_COMPONENTS)
          message(ERROR "Cannot build Boost.Python: your compiler is too old")
        endif()
      endif()
    endif()
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(BOOST_TOOLSET clang)
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(BOOST_TOOLSET intel-linux)
  endif()
  
  set(BOOST_WITH ${BOOST_COMPONENTS})
  string(REGEX REPLACE "[.]" "_" Boost_VERSION_NODOTS ${Boost_VERSION})
  if(BOOST_WITH)
    string(REPLACE ";" "," TMP "${BOOST_WITH}")
    set(BOOST_WITH "--with-libraries=${TMP}")
    list(APPEND BOOST_CONFIG_OPTS "${BOOST_WITH}")
  endif(BOOST_WITH)
  string(REPLACE ";" " " TMP "${BOOST_CONFIG_OPTS}")
  set(BOOST_CONFIG_OPTS "${TMP}")
  message(STATUS "BOOST_CONFIG_OPTS: ${BOOST_CONFIG_OPTS}")
  
  set(BOOST_URI http://sourceforge.net/projects/boost/files/boost/${Boost_VERSION}/boost_${Boost_VERSION_NODOTS}.tar.bz2/download)

  include(ExternalProject)
  ExternalProject_Add(external_boost
    URL ${BOOST_URI}
    DOWNLOAD_DIR ${DOWNLOAD_DIR}
    PREFIX ${BUILD_DIR}/boost
    CONFIGURE_COMMAND sh ./bootstrap.sh 
    --prefix=${BIN_DIR}/boost 
    --without-icu 
    --with-toolset=${BOOST_TOOLSET} ${BOOST_CONFIG_OPTS} && 
    echo "using mpi : ${MPI_CXX_COMPILER} $<SEMICOLON>" >> ./project-config.jam
    BUILD_COMMAND ./b2 -q ${BOOST_OPTS} stage
    INSTALL_COMMAND ./b2 -q ${BOOST_OPTS} install
    BUILD_IN_SOURCE 1
    URL_HASH MD5=${BOOST_MD5}
    )

  # There does not appear to be anything in Boost that refers to the
  # install path. So, things can just be copied. 
  
  install(
    DIRECTORY ${BIN_DIR}/boost/include/boost
    USE_SOURCE_PERMISSIONS
    DESTINATION include
    )

  install(
    DIRECTORY ${BIN_DIR}/boost/lib/
    DESTINATION lib
    USE_SOURCE_PERMISSIONS
    FILES_MATCHING PATTERN "*"
    )

  set(Boost_FOUND TRUE PARENT_SCOPE)
  set(Boost_INCLUDE_DIRS ${BIN_DIR}/boost/include PARENT_SCOPE)
  set(Boost_LIBRARY_DIRS ${BIN_DIR}/boost/lib PARENT_SCOPE)
  set(liblist "")
  message(STATUS "BOOST_COMPONENTS=${BOOST_COMPONENTS}")
  foreach(b IN LISTS BOOST_COMPONENTS)
    set(Boost_${b}_FOUND TRUE PARENT_SCOPE)
    list(APPEND liblist "${BIN_DIR}/boost/lib/libboost_${b}${LIB_SUFFIX}")
  endforeach()
  if(python IN_LIST BOOST_COMPONENTS)
    list(APPEND liblist ${PYTHON_LIBRARIES})
  endif()
  set(Boost_LIBRARIES "${liblist}" PARENT_SCOPE)
  
endfunction(BuildBoostExternalProject)