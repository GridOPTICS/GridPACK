// -------------------------------------------------------------
/**
 * @file   vector.cpp
 * @author William A. Perkins
 * @date   2013-05-10 12:12:22 d3g096
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

#include <iostream>
#include "gridpack/math/vector.hpp"
#include "gridpack/math/petsc/petsc_exception.hpp"
#include "gridpack/math/petsc/petsc_vector_implementation.hpp"
#include "gridpack/math/petsc/petsc_vector_extractor.hpp"

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
  p_vector_impl.reset(impl);
}

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
// vector clone
// -------------------------------------------------------------
Vector *
clone(const Vector& from)
{
  std::cerr << "In clone()" << std::endl;
  parallel::Communicator comm(from.communicator());
  int local_size(from.local_size());
  PETScVectorImplementation *pimpl =
    new PETScVectorImplementation(comm, local_size);

  const Vec* from_vec;
  {
    PETScConstVectorExtractor vext;
    from.accept(vext);
    from_vec = vext.vector();
  }

  Vec *to_vec;
  {
    PETScVectorExtractor vext;
    pimpl->accept(vext);
    to_vec = vext.vector();
  }

  PetscErrorCode ierr;

  try {
    ierr = VecCopy((*from_vec), *to_vec); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }

  Vector *result(new Vector(pimpl));
  
  return result;
}

} // namespace math
} // namespace gridpack
