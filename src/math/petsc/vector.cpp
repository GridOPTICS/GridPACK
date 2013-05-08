// -------------------------------------------------------------
/**
 * @file   vector.cpp
 * @author William A. Perkins
 * @date   2013-05-07 14:45:58 d3g096
 * 
 * @brief  PETSc implementation of gridpack::math::Vector
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
#include "gridpack/math/petsc/petsc_vector_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class Vector
// -------------------------------------------------------------

// -------------------------------------------------------------
// Vector:: constructors / destructor
// -------------------------------------------------------------
Vector::Vector(const parallel::Communicator& comm, const int& local_length)
  : parallel::Distributed(comm), utility::Uncopyable()
{
  PETScVectorImplementation *impl = 
    new PETScVectorImplementation(this->communicator(), local_length);
  vector_impl_.reset(impl);
}

Vector::~Vector(void)
{
  // empty
}

} // namespace math
} // namespace gridpack
