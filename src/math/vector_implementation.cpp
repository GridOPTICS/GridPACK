// -------------------------------------------------------------
/**
 * @file   vector_implementation.cpp
 * @author William A. Perkins
 * @date   2013-06-11 14:19:31 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "vector_implementation.hpp"

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

