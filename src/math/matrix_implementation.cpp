// -------------------------------------------------------------
/**
 * @file   matrix_implementation.cpp
 * @author William A. Perkins
 * @date   2013-06-11 14:16:58 d3g096
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
