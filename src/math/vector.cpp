// -------------------------------------------------------------
/**
 * @file   vector.cpp
 * @author William A. Perkins
 * @date   2013-05-15 13:38:27 d3g096
 * 
 * @brief  Part of Vector independent of specific implementation
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created May  7, 2013 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------

#include "gridpack/math/vector.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class Vector
// -------------------------------------------------------------

// -------------------------------------------------------------
// Vector:: constructors / destructor
// -------------------------------------------------------------
Vector::Vector(VectorImplementation *vimpl)
  : parallel::Distributed(vimpl->communicator()), utility::Uncopyable(),
    p_vector_impl(vimpl)
{
}

Vector::~Vector(void)
{
  // empty
}


} // namespace math
} // namespace gridpack
