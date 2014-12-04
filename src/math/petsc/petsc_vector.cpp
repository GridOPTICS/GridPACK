// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   vector.cpp
 * @author William A. Perkins
 * @date   2014-10-30 14:13:18 d3g096
 * 
 * @brief  PETSc-specific part of Vector
 * 
 * 
 */
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
  : parallel::WrappedDistributed(), utility::Uncopyable()
{
  PETScVectorImplementation<ComplexType> *impl = 
    new PETScVectorImplementation<ComplexType>(comm, local_length);
  p_vector_impl.reset(impl);
  p_setDistributed(impl);
}

// -------------------------------------------------------------
// Vector::add
// -------------------------------------------------------------
void
Vector::add(const Vector& x, const Vector::TheType& scale)
{
  this->p_checkCompatible(x);
  PetscErrorCode ierr(0);
  const Vec *xvec(PETScVector(x));
  Vec *yvec(PETScVector(*this));
  try {
    PetscScalar alpha(scale);

    // This call computes y = x + alpha*y. Where y is p_vector.  
    ierr = VecAXPY(*yvec, alpha, *xvec);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

void
Vector::add(const Vector::TheType& x)
{
  Vec *vec(PETScVector(*this));
  PetscErrorCode ierr(0);
  try {
    ierr = VecShift(*vec, x); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}


// -------------------------------------------------------------
// Vector::equate
// -------------------------------------------------------------
void
Vector::equate(const Vector& x)
{
  this->p_checkCompatible(x);
  PetscErrorCode ierr(0);
  Vec *yvec(PETScVector(*this));
  const Vec *xvec(PETScVector(x));
  try {
    ierr = VecCopy(*xvec, *yvec); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// Vector::elementMultiply
// -------------------------------------------------------------
void
Vector::elementMultiply(const Vector& x)
{
  Vec *vec(PETScVector(*this));
  const Vec *xvec(PETScVector(x));
  PetscErrorCode ierr(0);
  try {
    ierr = VecPointwiseMult(*vec, *vec, *xvec); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}  

// -------------------------------------------------------------
// Vector::elementDivide
// -------------------------------------------------------------
void
Vector::elementDivide(const Vector& x)
{
  Vec *vec(PETScVector(*this));
  const Vec *xvec(PETScVector(x));
  PetscErrorCode ierr(0);
  try {
    ierr = VecPointwiseDivide(*vec, *vec, *xvec); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}  


} // namespace math
} // namespace gridpack
