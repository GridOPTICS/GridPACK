#### Modified from http://www.openflipper.org/svnrepo/CoMISo/trunk/CoMISo/cmake/FindGUROBI.cmake


# - Try to find GUROBI
# Once done this will define
#  GUROBI_FOUND - System has Gurobi
#  GUROBI_INCLUDE_DIRS - The Gurobi include directories
#  GUROBI_LIBRARIES - The libraries needed to use Gurobi

if (GUROBI_INCLUDE_DIR)
  # in cache already
  set(GUROBI_FOUND TRUE)
  set(GUROBI_INCLUDE_DIRS ${GUROBI_INCLUDE_DIR} )
  set(GUROBI_LIBRARIES ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} )
else (GUROBI_INCLUDE_DIR)

find_path(GUROBI_INCLUDE_DIR 
          NAMES gurobi_c++.h
          PATHS ${GUROBI_HOME}/include
          )
message(STATUS "GUROBI Include Dir: ${GUROBI_INCLUDE_DIR}")

find_library( GUROBI_LIBRARY 
              NAMES gurobi
        gurobi56
        gurobi60        
              PATHS ${GUROBI_HOME}/lib 
              )
message(STATUS "GUROBI Library: ${GUROBI_LIBRARY}")

find_library( GUROBI_CXX_LIBRARY 
              NAMES gurobi_c++
              PATHS ${GUROBI_HOME}/lib 
              )
# use c++ headers as default
# handle the QUIETLY and REQUIRED arguments and set LIBGUROBI_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI  DEFAULT_MSG
                                  GUROBI_LIBRARY GUROBI_CXX_LIBRARY GUROBI_INCLUDE_DIR)


IF(GUROBI_FOUND)
  set(GUROBI_INCLUDE_DIRS ${GUROBI_INCLUDE_DIR})
  set(GUROBI_LIBRARIES ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY})
ENDIF(GUROBI_FOUND)

mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY GUROBI_CXX_LIBRARY)

endif(GUROBI_INCLUDE_DIR)
