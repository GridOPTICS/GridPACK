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
 * @date   2015-01-26 09:01:45 d3g096
 * 
 * @brief  PETSc-specific part of Vector
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include "vector.hpp"
#include "complex_operators.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_vector_implementation.hpp"
#include "petsc/petsc_vector_extractor.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class Vector
// -------------------------------------------------------------

// -------------------------------------------------------------
// VectorT:: constructors / destructor
// -------------------------------------------------------------
template <typename T, typename I>
VectorT<T, I>::VectorT(const parallel::Communicator& comm, const int& local_length)
  : parallel::WrappedDistributed(), utility::Uncopyable()
{
  PETScVectorImplementation<T, I> *impl = 
    new PETScVectorImplementation<T, I>(comm, local_length);
  p_vector_impl.reset(impl);
  p_setDistributed(impl);
}
template VectorT<ComplexType>::VectorT(const parallel::Communicator& comm, 
                                       const int& local_length);
template VectorT<RealType>::VectorT(const parallel::Communicator& comm, 
                                    const int& local_length);

// -------------------------------------------------------------
// VectorT::add
// -------------------------------------------------------------
template <typename T, typename I>
void
VectorT<T, I>::add(const VectorT<T, I>& x, const VectorT<T, I>::TheType& scale)
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

template void VectorT<ComplexType>::add(const VectorT<ComplexType>& x, 
                                        const typename VectorT<ComplexType>::TheType& scale);
template void VectorT<RealType>::add(const VectorT<RealType>& x, 
                                     const typename VectorT<RealType>::TheType& scale);

template <typename T, typename I>
void
VectorT<T, I>::add(const VectorT<T, I>::TheType& x)
{
  Vec *vec(PETScVector(*this));
  if (PETScVectorImplementation<T, I>::useLibrary::value) {
    PetscErrorCode ierr(0);
    try {
      PetscScalar tmpx =
        gridpack::math::equate<PetscScalar, TheType>(x);
      ierr = VecShift(*vec, x); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  } else {
    BOOST_ASSERT(false);
  }
}

template void VectorT<ComplexType, int>::add(const typename VectorT<ComplexType>::TheType& x);
template void VectorT<double, int>::add(const typename VectorT<RealType>::TheType& x);

// -------------------------------------------------------------
// VectorT::equate
// -------------------------------------------------------------
/** 
 * Make this vector equal to another.  This works regardless of the
 * underlying PETSc type. 
 * 
 * @param x vector to copy
 * 
 * @return 
 */
template <typename T, typename I>
void
VectorT<T, I>::equate(const VectorT<T, I>& x)
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

template void VectorT<ComplexType, int>::equate(const VectorT<ComplexType>& x);
template void VectorT<double, int>::equate(const VectorT<RealType>& x);

// -------------------------------------------------------------
// VectorT::elementMultiply
// -------------------------------------------------------------
template <typename T, typename I>
void
VectorT<T, I>::elementMultiply(const VectorT<T, I>& x)
{
  this->p_checkCompatible(x);
  Vec *vec(PETScVector(*this));
  const Vec *xvec(PETScVector(x));
  PetscErrorCode ierr(0);
  try {
    if (PETScVectorImplementation<T, I>::useLibrary::value) {
      ierr = VecPointwiseMult(*vec, *vec, *xvec); CHKERRXX(ierr);
    } else {
      BOOST_ASSERT(false);
    }
  } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
  }
}  

template void VectorT<ComplexType, int>::elementMultiply(const VectorT<ComplexType>& x);
template void VectorT<double, int>::elementMultiply(const VectorT<RealType>& x);

// -------------------------------------------------------------
// VectorT::elementDivide
// -------------------------------------------------------------
template <typename T, typename I>
void
VectorT<T, I>::elementDivide(const VectorT<T, I>& x)
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

template void VectorT<ComplexType, int>::elementDivide(const VectorT<ComplexType> &x);
template void VectorT<double, int>::elementDivide(const VectorT<RealType> &x);


} // namespace math
} // namespace gridpack
