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
 * @date   2015-06-12 14:01:18 d3g096
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
// applyBinaryOperator
// -------------------------------------------------------------
template <typename T>
static void
applyBinaryOperator(Vec *v1, const Vec *v2, 
                    base_binary_function<T>& op)
{
  PetscErrorCode ierr;

  PetscInt n1, n2;
  ierr = VecGetLocalSize(*v1, &n1); CHKERRXX(ierr);
  ierr = VecGetLocalSize(*v2, &n2); CHKERRXX(ierr);
  BOOST_ASSERT(n1 == n2);

  PetscScalar *x1;
  const PetscScalar *x2;
  ierr = VecGetArray(*v1, &x1);  CHKERRXX(ierr);
  ierr = VecGetArrayRead(*v2, &x2); CHKERRXX(ierr);
  binary_operation<T, PetscScalar>(n1, x1, x2, op);
  ierr = VecRestoreArrayRead(*v2, &x2); CHKERRXX(ierr);
  ierr = VecRestoreArray(*v1, &x1);  CHKERRXX(ierr);
  ierr = VecAssemblyBegin(*v1); CHKERRXX(ierr);
  ierr = VecAssemblyEnd(*v1); 
}

// -------------------------------------------------------------
// applyBinaryOperator
// -------------------------------------------------------------
template <typename T>
static void
applyUnaryOperator(Vec *v, base_unary_function<T>& op)
{
  PetscErrorCode ierr;

  PetscInt n;
  ierr = VecGetLocalSize(*v, &n); CHKERRXX(ierr);
  PetscScalar *x;
  ierr = VecGetArray(*v, &x);  CHKERRXX(ierr);
  unary_operation<T, PetscScalar>(n, x, op);
  ierr = VecRestoreArray(*v, &x);  CHKERRXX(ierr);
  ierr = VecAssemblyBegin(*v); CHKERRXX(ierr);
  ierr = VecAssemblyEnd(*v); 
}

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
/** 
 * This should work fine regardless of what PetscScalar is.
 * 
 * @param T 
 * @param x 
 * @param T 
 * @param scale 
 * 
 * @return 
 */
template <typename T, typename I>
void
VectorT<T, I>::add(const VectorT<T, I>& x, const VectorT<T, I>::TheType& scale)
{
  this->p_checkCompatible(x);
  
  PetscErrorCode ierr(0);
  const Vec *xvec(PETScVector(x));
  Vec *yvec(PETScVector(*this));
  try {
    if (PETScVectorImplementation<T, I>::useLibrary) {
      PetscScalar alpha =
        gridpack::math::equate<PetscScalar, TheType>(scale);
      // This call computes y = x + alpha*y. Where y is p_vector.  
      ierr = VecAXPY(*yvec, alpha, *xvec);
    } else {
      scaleAdd2<T> op(scale);
      applyBinaryOperator<T>(yvec, xvec, op);
    }
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

template void VectorT<ComplexType>::add(const VectorT<ComplexType>& x, 
                                        const VectorT<ComplexType>::TheType& scale);
template void VectorT<RealType>::add(const VectorT<RealType>& x, 
                                     const VectorT<RealType>::TheType& scale);

template <typename T, typename I>
void
VectorT<T, I>::add(const VectorT<T, I>::TheType& x)
{
  if (PETScVectorImplementation<T, I>::useLibrary) {
    Vec *vec(PETScVector(*this));
    PetscErrorCode ierr(0);
    try {
      if (PETScVectorImplementation<T, I>::useLibrary) {
        PetscScalar tmpx =
          gridpack::math::equate<PetscScalar, TheType>(x);
        ierr = VecShift(*vec, tmpx); CHKERRXX(ierr);
      } else {
        addvalue<TheType> op(x);
        applyUnaryOperator(vec, op);
      }        
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  } else {
    gridpack::math::addvalue<TheType> op(x);
    this->p_vector_impl->applyOperation(op);
  }
}

template void VectorT<ComplexType, int>::add(const VectorT<ComplexType>::TheType& x);
template void VectorT<double, int>::add(const VectorT<RealType>::TheType& x);

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
    if (PETScVectorImplementation<T, I>::useLibrary) {
      ierr = VecPointwiseMult(*vec, *vec, *xvec); CHKERRXX(ierr);
    } else {
      multiplyvalue2<T> op;
      applyBinaryOperator<T>(vec, xvec, op);
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
  this->p_checkCompatible(x);
  Vec *vec(PETScVector(*this));
  const Vec *xvec(PETScVector(x));
  PetscErrorCode ierr(0);
  try {
    if (PETScVectorImplementation<T, I>::useLibrary) {
      ierr = VecPointwiseDivide(*vec, *vec, *xvec); CHKERRXX(ierr);
    } else {
      dividevalue2<T> op;
      applyBinaryOperator<T>(vec, xvec, op);
    }
  } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
  }
}  

template void VectorT<ComplexType, int>::elementDivide(const VectorT<ComplexType> &x);
template void VectorT<double, int>::elementDivide(const VectorT<RealType> &x);


} // namespace math
} // namespace gridpack
