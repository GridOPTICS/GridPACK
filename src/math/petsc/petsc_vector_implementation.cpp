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
 * @date   2013-10-09 13:25:41 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


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
                                                     const int& local_length)
  : VectorImplementation(comm), p_min_index(-1), p_max_index(-1)
{
  PetscErrorCode ierr;
  try {

    // If any ownership arguments are specifed, *all* ownership arguments
    // need to be specified.

    PetscInt llen(local_length), glen(PETSC_DETERMINE);
    ierr = PetscSplitOwnership(comm, &llen, &glen); CHKERRXX(ierr);


    ierr = VecCreate(comm,&p_vector); CHKERRXX(ierr);
    ierr = VecSetSizes(p_vector, llen, glen); CHKERRXX(ierr);
    if (comm.size() > 1) {
      ierr = VecSetType(p_vector, VECMPI);  CHKERRXX(ierr);
    } else {
      ierr = VecSetType(p_vector, VECSEQ);  CHKERRXX(ierr);
    }

    // set and gets only work for values on this processor
    ierr = VecSetOption(p_vector, VEC_IGNORE_OFF_PROC_ENTRIES, PetscBool(true)); CHKERRXX(ierr);

    // get and save the ownership index range
    PetscInt lo, hi;
    ierr = VecGetOwnershipRange(p_vector, &lo, &hi); CHKERRXX(ierr);
    p_min_index = lo;
    p_max_index = hi;

  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

PETScVectorImplementation::~PETScVectorImplementation(void)
{
  // Bad things happen (e.g. race condition on RHEL5) if one tries
  // to destroy a PETSc thing after PETSc is finalized.
  PetscErrorCode ierr;
  
  try  {
    PetscBool ok;
    ierr = PetscInitialized(&ok);
    if (ok) {
      ierr = VecDestroy(&p_vector);
    }
  } catch (...) {
    // just eat it
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_size
// -------------------------------------------------------------
int 
PETScVectorImplementation::p_size(void) const
{
  PetscErrorCode ierr;
  try {
    PetscInt gsize;
    ierr = VecGetSize(this->p_vector, &gsize); CHKERRXX(ierr);
    return static_cast<int>(gsize);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_localSize
// -------------------------------------------------------------
int 
PETScVectorImplementation::p_localSize(void) const
{
  PetscErrorCode ierr;
  try {
    PetscInt gsize;
    ierr = VecGetLocalSize(this->p_vector, &gsize); CHKERRXX(ierr);
    return static_cast<int>(gsize);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_localIndexRange
// -------------------------------------------------------------
void
PETScVectorImplementation::p_localIndexRange(int& lo, int& hi) const
{
  lo = p_min_index;
  hi = p_max_index;
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
PETScVectorImplementation::p_setElement(const int& i, const ComplexType& x)
{
  PetscErrorCode ierr;
  try {
    ierr = VecSetValue(p_vector, i, x, INSERT_VALUES); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_setElements
// -------------------------------------------------------------
void
PETScVectorImplementation::p_setElements(const int& n, const int *i, const ComplexType *x)
{
  PetscErrorCode ierr;
  try {
    ierr = VecSetValues(p_vector, n, i, x, INSERT_VALUES); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_addElement
// -------------------------------------------------------------
void
PETScVectorImplementation::p_addElement(const int& i, const ComplexType& x)
{
  PetscErrorCode ierr;
  try {
    ierr = VecSetValue(p_vector, i, x, ADD_VALUES); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_addElements
// -------------------------------------------------------------
void
PETScVectorImplementation::p_addElements(const int& n, const int *i, const ComplexType *x)
{
  PetscErrorCode ierr;
  try {
    ierr = VecSetValues(p_vector, n, i, x, ADD_VALUES); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_getElement
// -------------------------------------------------------------
void
PETScVectorImplementation::p_getElement(const int& i, ComplexType& x) const
{
  this->p_getElements(1, &i, &x);
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_get_elements
// -------------------------------------------------------------
void
PETScVectorImplementation::p_getElements(const int& n, const int *i, ComplexType *x) const
{
  // FIXME: Cannot get off process elements
  PetscErrorCode ierr;
  try {
    ierr = VecGetValues(p_vector, n, i, x); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_get_all_elements
// -------------------------------------------------------------
void
PETScVectorImplementation::p_getAllElements(ComplexType *x) const
{
  PetscErrorCode ierr(0);
  try {
    VecScatter scatter;
    Vec all;
    int n(this->size());
    ierr = VecScatterCreateToAll(p_vector, &scatter, &all); CHKERRXX(ierr);
    ierr = VecScatterBegin(scatter, p_vector, all, INSERT_VALUES, SCATTER_FORWARD); CHKERRXX(ierr);
    ierr = VecScatterEnd(scatter, p_vector, all, INSERT_VALUES, SCATTER_FORWARD); CHKERRXX(ierr);
    const PetscScalar *tmp;
    ierr = VecGetArrayRead(all, &tmp); CHKERRXX(ierr);
    std::copy(tmp, tmp + n, &x[0]);
    ierr = VecRestoreArrayRead(all, &tmp); CHKERRXX(ierr);
    ierr = VecScatterDestroy(&scatter); CHKERRXX(ierr);
    ierr = VecDestroy(&all); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_zero
// -------------------------------------------------------------
void
PETScVectorImplementation::p_zero(void)
{
  ComplexType v(0.0, 0.0);
  this->fill(v);
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_fill
// -------------------------------------------------------------
void
PETScVectorImplementation::p_fill(const ComplexType& v)
{
  PetscErrorCode ierr(0);
  try {
    PetscScalar pv(v);
    ierr = VecSet(this->p_vector, pv); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}  

// -------------------------------------------------------------
// PETScVectorImplementation::p_norm1
// -------------------------------------------------------------
ComplexType
PETScVectorImplementation::p_norm1(void) const
{
  ComplexType result;
  PetscErrorCode ierr(0);
  try {
    PetscReal v;
    ierr = VecNorm(this->p_vector, NORM_1, &v); CHKERRXX(ierr);
    result = v;
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  return result;
}


// -------------------------------------------------------------
// PETScVectorImplementation::p_norm2
// -------------------------------------------------------------
ComplexType
PETScVectorImplementation::p_norm2(void) const
{
  ComplexType result;
  PetscErrorCode ierr(0);
  try {
    PetscReal v;
    ierr = VecNorm(this->p_vector, NORM_2, &v); CHKERRXX(ierr);
    result = v;
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  return result;
}


// -------------------------------------------------------------
// PETScVectorImplementation::p_ready
// -------------------------------------------------------------
void
PETScVectorImplementation::p_ready(void)
{
  PetscErrorCode ierr;
  try {
    ierr = VecAssemblyBegin(p_vector); CHKERRXX(ierr);
    ierr = VecAssemblyEnd(p_vector); 
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }

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
  int local_size(this->localSize());
  
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
    ierr = VecCopy(this->p_vector, *to_vec); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  return result;
}


} // namespace math
} // namespace gridpack

