// -------------------------------------------------------------
/**
 * @file   vector.cpp
 * @author William A. Perkins
 * @date   2013-06-11 14:12:28 d3g096
 * 
 * @brief  PETSc-specific part of Vector
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
#include "vector.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_vector_implementation.hpp"
#include "petsc/petsc_vector_extractor.hpp"

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

// -------------------------------------------------------------
// Vector::scale
// -------------------------------------------------------------
void 
Vector::scale(const ComplexType& x)
{
  Vec *vec(PETScVector(*this));
  PetscErrorCode ierr(0);
  try {
    ierr = VecScale(*vec, x); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// Vector::add
// -------------------------------------------------------------
void
Vector::add(const Vector& x)
{
  this->p_check_compatible(x);
  PetscErrorCode ierr(0);
  const Vec *xvec(PETScVector(x));
  Vec *yvec(PETScVector(*this));
  try {
    PetscScalar alpha(1.0);

    // This call computes y = x + alpha*y. Where y is p_vector.  
    ierr = VecAYPX(*yvec, alpha, *xvec);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

void
Vector::add(const ComplexType& x)
{
  Vec *vec(PETScVector(*this));
  PetscErrorCode ierr(0);
  try {
    ierr = VecShift(*vec, x); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}





// -------------------------------------------------------------
// Vector::equate
// -------------------------------------------------------------
void
Vector::equate(const Vector& x)
{
  this->p_check_compatible(x);
  PetscErrorCode ierr(0);
  Vec *yvec(PETScVector(*this));
  const Vec *xvec(PETScVector(x));
  try {
    ierr = VecCopy(*xvec, *yvec); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// Vector::reciprocal
// -------------------------------------------------------------
void
Vector::reciprocal(void)
{
  Vec *vec(PETScVector(*this));
  PetscErrorCode ierr(0);
  try {
    ierr = VecReciprocal(*vec); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}




} // namespace math
} // namespace gridpack
