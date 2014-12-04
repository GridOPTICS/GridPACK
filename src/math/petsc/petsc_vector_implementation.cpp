// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_vector_implementation.cpp
 * @author William A. Perkins
 * @date   2014-10-28 14:26:47 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <petscsys.h>
#include "implementation_visitor.hpp"
#include "petsc/petsc_vector_implementation.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_vector_extractor.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScVectorImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScVectorImplementation:: constructors / destructor
// -------------------------------------------------------------
PETScVectorImplementation::PETScVectorImplementation(const parallel::Communicator& comm,
                                                     const PETScVectorImplementation::IdxType& local_length)
  : VectorImplementation(comm), p_vwrap(comm, local_length)
{
}

PETScVectorImplementation::PETScVectorImplementation(Vec& pVec, const bool& copyVec)
  : VectorImplementation(PetscVectorWrapper::getCommunicator(pVec)), 
    p_vwrap(pVec, copyVec)
{
}


PETScVectorImplementation::~PETScVectorImplementation(void)
{
}

// -------------------------------------------------------------
// PETScVectorImplementation::getVector
// -------------------------------------------------------------
const Vec *
PETScVectorImplementation::getVector(void) const
{
  return p_vwrap.getVector();
}

Vec *
PETScVectorImplementation::getVector(void)
{
  return p_vwrap.getVector();
}


// -------------------------------------------------------------
// PETScVectorImplementation::p_size
// -------------------------------------------------------------
PETScVectorImplementation::IdxType 
PETScVectorImplementation::p_size(void) const
{
  return p_vwrap.size();
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_localSize
// -------------------------------------------------------------
PETScVectorImplementation::IdxType 
PETScVectorImplementation::p_localSize(void) const
{
  return p_vwrap.localSize();
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_localIndexRange
// -------------------------------------------------------------
void
PETScVectorImplementation::p_localIndexRange(PETScVectorImplementation::IdxType& lo, 
                                             PETScVectorImplementation::IdxType& hi) const
{
  PetscInt plo, phi;
  p_vwrap.localIndexRange(plo, phi);
  lo = plo;
  hi = phi;
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_setElement
// -------------------------------------------------------------
/** 
 * If you try to set an off processor value, it will be ignored
 * 
 * @param i global vector index
 * @param x value
 */
void
PETScVectorImplementation::p_setElement(const PETScVectorImplementation::IdxType& i, 
                                        const PETScVectorImplementation::TheType& x)
{
  PetscErrorCode ierr;
  try {
    Vec *v = p_vwrap.getVector();
    ierr = VecSetValue(*v, i, x, INSERT_VALUES); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_setElements
// -------------------------------------------------------------
void
PETScVectorImplementation::p_setElements(const PETScVectorImplementation::IdxType& n, 
                                         const PETScVectorImplementation::IdxType *i, 
                                         const PETScVectorImplementation::TheType *x)
{
  PetscErrorCode ierr;
  try {
    Vec *v = p_vwrap.getVector();
    ierr = VecSetValues(*v, n, i, x, INSERT_VALUES); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_addElement
// -------------------------------------------------------------
void
PETScVectorImplementation::p_addElement(const PETScVectorImplementation::IdxType& i, 
                                        const PETScVectorImplementation::TheType& x)
{
  PetscErrorCode ierr;
  try {
    Vec *v = p_vwrap.getVector();
    ierr = VecSetValue(*v, i, x, ADD_VALUES); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_addElements
// -------------------------------------------------------------
void
PETScVectorImplementation::p_addElements(const PETScVectorImplementation::IdxType& n, 
                                         const PETScVectorImplementation::IdxType *i, 
                                         const PETScVectorImplementation::TheType *x)
{
  PetscErrorCode ierr;
  try {
   Vec *v = p_vwrap.getVector();
   ierr = VecSetValues(*v, n, i, x, ADD_VALUES); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_getElement
// -------------------------------------------------------------
void
PETScVectorImplementation::p_getElement(const PETScVectorImplementation::IdxType& i, 
                                        PETScVectorImplementation::TheType& x) const
{
  this->p_getElements(1, &i, &x);
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_get_elements
// -------------------------------------------------------------
void
PETScVectorImplementation::p_getElements(const PETScVectorImplementation::IdxType& n, 
                                         const PETScVectorImplementation::IdxType *i, 
                                         PETScVectorImplementation::TheType *x) const
{
  // FIXME: Cannot get off process elements
  PetscErrorCode ierr;
  try {
    const Vec *v = p_vwrap.getVector();
    ierr = VecGetValues(*v, n, i, x); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_get_all_elements
// -------------------------------------------------------------
void
PETScVectorImplementation::p_getAllElements(PETScVectorImplementation::TheType *x) const
{
  PetscErrorCode ierr(0);
  try {
    const Vec *v = p_vwrap.getVector();
    VecScatter scatter;
    Vec all;
    IdxType n(this->size());
    ierr = VecScatterCreateToAll(*v, &scatter, &all); CHKERRXX(ierr);
    ierr = VecScatterBegin(scatter, *v, all, INSERT_VALUES, SCATTER_FORWARD); CHKERRXX(ierr);
    ierr = VecScatterEnd(scatter, *v, all, INSERT_VALUES, SCATTER_FORWARD); CHKERRXX(ierr);
    const PetscScalar *tmp;
    ierr = VecGetArrayRead(all, &tmp); CHKERRXX(ierr);
    std::copy(tmp, tmp + n, &x[0]);
    ierr = VecRestoreArrayRead(all, &tmp); CHKERRXX(ierr);
    ierr = VecScatterDestroy(&scatter); CHKERRXX(ierr);
    ierr = VecDestroy(&all); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_zero
// -------------------------------------------------------------
void
PETScVectorImplementation::p_zero(void)
{
  p_vwrap.zero();
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_fill
// -------------------------------------------------------------
void
PETScVectorImplementation::p_fill(const PETScVectorImplementation::TheType& v)
{
  PetscErrorCode ierr(0);
  try {
    PetscScalar pv(v);
    Vec *v = p_vwrap.getVector();
    ierr = VecSet(*v, pv); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}  

// -------------------------------------------------------------
// PETScVectorImplementation::p_norm1
// -------------------------------------------------------------
double
PETScVectorImplementation::p_norm1(void) const
{
  return p_vwrap.norm1();
}


// -------------------------------------------------------------
// PETScVectorImplementation::p_norm2
// -------------------------------------------------------------
double
PETScVectorImplementation::p_norm2(void) const
{
  return p_vwrap.norm2();
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_normInfinity
// -------------------------------------------------------------
double
PETScVectorImplementation::p_normInfinity(void) const
{
  return p_vwrap.normInfinity();
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_abs
// -------------------------------------------------------------
void
PETScVectorImplementation::p_abs(void)
{
  p_vwrap.abs();
}
  

// -------------------------------------------------------------
// PETScVectorImplementation::p_conjugate
// -------------------------------------------------------------
void
PETScVectorImplementation::p_conjugate(void)
{
  p_vwrap.conjugate();
}


// -------------------------------------------------------------
// PETScVectorImplementation::p_exp
// -------------------------------------------------------------
void
PETScVectorImplementation::p_exp(void)
{
  p_vwrap.exp();
}



// -------------------------------------------------------------
// PETScVectorImplementation::p_ready
// -------------------------------------------------------------
void
PETScVectorImplementation::p_ready(void)
{
  p_vwrap.ready();
}


// -------------------------------------------------------------
// PETScVectorImplementation::p_accept
// -------------------------------------------------------------
void 
PETScVectorImplementation::p_accept(ImplementationVisitor& visitor)
{
  visitor.visit(*this);
}

void 
PETScVectorImplementation::p_accept(ConstImplementationVisitor& visitor) const
{
  visitor.visit(*this);
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_clone
// -------------------------------------------------------------
VectorImplementation *
 PETScVectorImplementation::p_clone(void) const
{
  parallel::Communicator comm(this->communicator());
  IdxType local_size(this->localSize());
  
  PETScVectorImplementation *result = 
    new PETScVectorImplementation(comm, local_size);
  PetscErrorCode ierr;

  Vec *to_vec;
  {
    PETScVectorExtractor vext;
    result->accept(vext);
    to_vec = vext.vector();
  }


  try {
    const Vec *v = p_vwrap.getVector();
    ierr = VecCopy(*v, *to_vec); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
  return result;
}


} // namespace math
} // namespace gridpack

