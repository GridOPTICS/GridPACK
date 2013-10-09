// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   matrix_implementation.cpp
 * @author William A. Perkins
 * @date   2013-10-09 12:24:04 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------

#include "matrix_implementation.hpp"

namespace gridpack {
namespace math {


// -------------------------------------------------------------
//  class MatrixImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// MatrixImplementation:: constructors / destructor
// -------------------------------------------------------------
MatrixImplementation::MatrixImplementation(const parallel::Communicator& comm)
  : utility::Uncopyable(), parallel::Distributed(comm)
{
  
}

MatrixImplementation::~MatrixImplementation(void)
{
}

} // namespace math
} // namespace gridpack
