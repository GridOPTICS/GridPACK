// -------------------------------------------------------------
/**
 * @file   vector_implementation.cpp
 * @author William A. Perkins
 * @date   2013-05-08 15:19:44 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/math/vector_implementation.hpp"

namespace gridpack {
namespace math {


// -------------------------------------------------------------
//  class VectorImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// VectorImplementation:: constructors / destructor
// -------------------------------------------------------------
VectorImplementation::VectorImplementation(const parallel::Communicator& comm)
  : utility::Uncopyable(), parallel::Distributed(comm)
{
  
}

VectorImplementation::~VectorImplementation(void)
{
}

} // namespace math
} // namespace gridpack

