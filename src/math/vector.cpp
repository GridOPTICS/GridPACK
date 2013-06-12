// -------------------------------------------------------------
/**
 * @file   vector.cpp
 * @author William A. Perkins
 * @date   2013-06-11 14:11:51 d3g096
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

#include "gridpack/utilities/exception.hpp"
#include "vector.hpp"

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

// -------------------------------------------------------------
// Vector::p_check_compatible
// -------------------------------------------------------------
void
Vector::p_check_compatible(const Vector& x) const
{
  if (this->communicator() != x.communicator()) {
    throw gridpack::Exception("incompatible: communicators do not match");
  }

  if (this->size() != x.size()) {
    throw gridpack::Exception("incompatible: sizes do not match");
  }
}

// -------------------------------------------------------------
// add
// -------------------------------------------------------------
Vector *
add(const Vector& A, const Vector& B)
{
  Vector *result(A.clone());
  result->add(B);
  return result;
}

void
add(const Vector& A, const Vector& B, Vector& result)
{
  result.equate(A);
  result.add(B);
}

} // namespace math
} // namespace gridpack
