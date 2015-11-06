# This module finds cplex.
#
# User can give CPLEX_ROOT_DIR as a hint stored in the cmake cache.
#
# It sets the following variables:
# CPLEX_FOUND - Set to false, or undefined, if cplex isn't found.
# CPLEX_INCLUDE_DIRS - include directory
# CPLEX_LIBRARIES - library files

FIND_PATH(CPLEX_INCLUDE_DIR ilcplex/cplex.h
  HINTS ${CPLEX_ROOT_DIR}/cplex/include
        ${CPLEX_ROOT_DIR}/include
  PATHS ENV C_INCLUDE_PATH
        ENV C_PLUS_INCLUDE_PATH
        ENV INCLUDE_PATH
)

FIND_PATH(CPLEX_CONCERT_INCLUDE_DIR
  ilconcert/iloenv.h
  HINTS ${CPLEX_ROOT_DIR}/concert/include
        ${CPLEX_ROOT_DIR}/include
  PATHS ENV C_INCLUDE_PATH
        ENV C_PLUS_INCLUDE_PATH
        ENV INCLUDE_PATH
)

FIND_LIBRARY(CPLEX_LIBRARY
  NAMES cplex${CPLEX_WIN_VERSION} cplex
  HINTS ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_debian4.0_4.1/static_pic
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_sles10_4.1/static_pic
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_linux/static_pic
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_osx/static_pic
  PATHS ENV LIBRARY_PATH
        ENV LD_LIBRARY_PATH
)
message(STATUS "CPLEX Library: ${CPLEX_LIBRARY}")

FIND_LIBRARY(CPLEX_ILOCPLEX_LIBRARY
  ilocplex
  HINTS ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_debian4.0_4.1/static_pic #unix
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_sles10_4.1/static_pic #unix
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_linux/static_pic
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_osx/static_pic
  PATHS ENV LIBRARY_PATH
        ENV LD_LIBRARY_PATH
)
message(STATUS "ILOCPLEX Library: ${CPLEX_ILOCPLEX_LIBRARY}")

FIND_LIBRARY(CPLEX_CONCERT_LIBRARY
  concert
  HINTS ${CPLEX_ROOT_DIR}/concert/lib/x86-64_debian4.0_4.1/static_pic #unix
        ${CPLEX_ROOT_DIR}/concert/lib/x86-64_sles10_4.1/static_pic #unix
        ${CPLEX_ROOT_DIR}/concert/lib/x86-64_linux/static_pic 
        ${CPLEX_ROOT_DIR}/concert/lib/x86-64_osx/static_pic
  PATHS ENV LIBRARY_PATH
        ENV LD_LIBRARY_PATH
)
message(STATUS "CONCERT Library: ${CPLEX_CONCERT_LIBRARY}")

FIND_PATH(CPLEX_BIN_DIR
  cplex
  HINTS ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_sles10_4.1 #unix
  ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_debian4.0_4.1 #unix
  ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_linux #unix
  ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_osx
  ENV LIBRARY_PATH
  ENV LD_LIBRARY_PATH
)
message(STATUS "CPLEX Bin Dir: ${CPLEX_BIN_DIR}")
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPLEX DEFAULT_MSG
  CPLEX_LIBRARY CPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_INCLUDE_DIR)

IF(CPLEX_FOUND)
  SET(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIR})
  SET(CPLEX_LIBRARIES ${CPLEX_CONCERT_LIBRARY} ${CPLEX_ILOCPLEX_LIBRARY} ${CPLEX_LIBRARY} )
  IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    SET(CPLEX_LIBRARIES "${CPLEX_LIBRARIES};m;pthread")
  ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
ENDIF(CPLEX_FOUND)

MARK_AS_ADVANCED(CPLEX_LIBRARY CPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY CPLEX_CONCERT_INCLUDE_DIR CPLEX_CONCERT_LIBRARY)
