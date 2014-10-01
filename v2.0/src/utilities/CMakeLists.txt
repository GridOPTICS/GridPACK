#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be
#     found
#     in the LICENSE file in the top level directory of this distribution.
#
# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created June 11, 2013 by William A. Perkins
# Last Change: 2013-06-11 15:00:23 d3g096
# -------------------------------------------------------------

# -------------------------------------------------------------
# gridpack header installation
# -------------------------------------------------------------
install(FILES 
  complex.hpp
  exception.hpp
  uncopyable.hpp
  named.hpp
  doxygen_dummy.hpp
  indexed.hpp
  sanity_interface.hpp
  DESTINATION include/gridpack/utilities
)

