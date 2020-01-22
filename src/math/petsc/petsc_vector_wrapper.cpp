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
 * @date   2019-09-13 13:56:23 d3g096
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
#include "petsc/petsc_vector_wrapper.hpp"
#include "implementation_visitor.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PetscVectorWrapper
// -------------------------------------------------------------

// -------------------------------------------------------------
// static PetscVectorWrapper::getCommunicator
// -------------------------------------------------------------
parallel::Communicator
PetscVectorWrapper::getCommunicator(const Vec& v)
{
  MPI_Comm comm(PetscObjectComm((PetscObject)v));
  parallel::Communicator result(comm);
  return result;
}



// -------------------------------------------------------------
// PetscVectorWrapper:: constructors / destructor
// -------------------------------------------------------------
PetscVectorWrapper::PetscVectorWrapper(const parallel::Communicator& comm,
                                       const PetscInt& local_length)
  : p_minIndex(-1), p_maxIndex(-1), p_vectorWrapped(false)
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
    p_minIndex = lo;
    p_maxIndex = hi;

  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
  
}

PetscVectorWrapper::PetscVectorWrapper(const PetscVectorWrapper& old)
  : p_minIndex(-1), p_maxIndex(-1), p_vectorWrapped(false)
{
  PetscErrorCode ierr;
  try {
    ierr = VecCopy(old.p_vector, p_vector); CHKERRXX(ierr);

    // get and save the ownership index range
    PetscInt lo, hi;
    ierr = VecGetOwnershipRange(p_vector, &lo, &hi); CHKERRXX(ierr);
    p_minIndex = lo;
    p_maxIndex = hi;

  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

PetscVectorWrapper::PetscVectorWrapper(Vec& pVec, const bool& copyVec)
  :p_minIndex(-1), p_maxIndex(-1), p_vectorWrapped(false)
{
  PetscErrorCode ierr;
  try {

    if (copyVec) {
      ierr = VecCopy(pVec, p_vector); CHKERRXX(ierr);
    } else {
      p_vector = pVec;
      p_vectorWrapped = true;
    }

    // get and save the ownership index range
    PetscInt lo, hi;
    ierr = VecGetOwnershipRange(p_vector, &lo, &hi); CHKERRXX(ierr);
    p_minIndex = lo;
    p_maxIndex = hi;

  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

PetscVectorWrapper::~PetscVectorWrapper(void)
{
  // Bad things happen (e.g. race condition on RHEL5) if one tries
  // to destroy a PETSc thing after PETSc is finalized.
  PetscErrorCode ierr(0);
  
  if (!p_vectorWrapped) {
    try  {
      PetscBool ok;
      ierr = PetscInitialized(&ok); CHKERRXX(ierr);
      if (ok) {
        ierr = VecDestroy(&p_vector); CHKERRXX(ierr);
      }
    } catch (...) {
      // just eat it
    }
  }
}

// -------------------------------------------------------------
// PetscVectorWrapper::size
// -------------------------------------------------------------
PetscInt 
PetscVectorWrapper::size(void) const
{
  PetscErrorCode ierr;
  try {
    PetscInt gsize;
    ierr = VecGetSize(this->p_vector, &gsize); CHKERRXX(ierr);
    return static_cast<PetscInt>(gsize);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscVectorWrapper::localSize
// -------------------------------------------------------------
PetscInt 
PetscVectorWrapper::localSize(void) const
{
  PetscErrorCode ierr;
  try {
    PetscInt gsize;
    ierr = VecGetLocalSize(this->p_vector, &gsize); CHKERRXX(ierr);
    return static_cast<PetscInt>(gsize);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscVectorWrapper::localIndexRange
// -------------------------------------------------------------
void
PetscVectorWrapper::localIndexRange(PetscInt& lo, 
                                      PetscInt& hi) const
{
  lo = p_minIndex;
  hi = p_maxIndex;
}

// -------------------------------------------------------------
// PetscVectorWrapper::getAllElements
// -------------------------------------------------------------
void
PetscVectorWrapper::getAllElements(PetscScalar *x) const
{
  PetscErrorCode ierr(0);
  try {
    const Vec *v = getVector();
    VecScatter scatter;
    Vec all;
    PetscInt n(this->size());
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
// PetscVectorWrapper::zero
// -------------------------------------------------------------
void
PetscVectorWrapper::zero(void)
{
  PetscErrorCode ierr;
  try {
    ierr = VecZeroEntries(this->p_vector); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}


// -------------------------------------------------------------
// PetscVectorWrapper::p_norm
// -------------------------------------------------------------
double
PetscVectorWrapper::p_norm(const NormType& t) const
{
  double result;
  PetscErrorCode ierr(0);
  try {
    PetscReal v;
    ierr = VecNorm(this->p_vector, t, &v); CHKERRXX(ierr);
    result = v;
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
  return result;
}


// -------------------------------------------------------------
// PetscVectorWrapper::norm1
// -------------------------------------------------------------
double
PetscVectorWrapper::norm1(void) const
{
  return p_norm(NORM_1);
}


// -------------------------------------------------------------
// PetscVectorWrapper::norm2
// -------------------------------------------------------------
double
PetscVectorWrapper::norm2(void) const
{
  return p_norm(NORM_2);
}

// -------------------------------------------------------------
// PetscVectorWrapper::normInfinity
// -------------------------------------------------------------
double
PetscVectorWrapper::normInfinity(void) const
{
  return p_norm(NORM_INFINITY);
}

// -------------------------------------------------------------
// PetscVectorWrapper::abs
// -------------------------------------------------------------
void
PetscVectorWrapper::abs(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = VecAbs(this->p_vector); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}  

// -------------------------------------------------------------
// PetscVectorWrapper::conjugate
// -------------------------------------------------------------
void
PetscVectorWrapper::conjugate(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = VecConjugate(this->p_vector); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}  

// -------------------------------------------------------------
// PetscVectorWrapper::exp
// -------------------------------------------------------------
void
PetscVectorWrapper::exp(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = VecExp(this->p_vector); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}  

// -------------------------------------------------------------
// PetscVectorWrapper::reciprocal
// -------------------------------------------------------------
void
PetscVectorWrapper::reciprocal(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = VecReciprocal(p_vector); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscVectorWrapper::ready
// -------------------------------------------------------------
void
PetscVectorWrapper::ready(void)
{
  PetscErrorCode ierr;
  try {
    ierr = VecAssemblyBegin(p_vector); CHKERRXX(ierr);
    ierr = VecAssemblyEnd(p_vector); 
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }

}

// -------------------------------------------------------------
// petsc_print_vector
// -------------------------------------------------------------
static void
petsc_print_vector(const Vec vec, const char* filename, PetscViewerFormat format)
{
  PetscErrorCode ierr;
  try {
    PetscViewer viewer;
    MPI_Comm comm = PetscObjectComm((PetscObject)vec);
    if (filename != NULL) {
      ierr = PetscViewerASCIIOpen(comm, filename, &viewer); ; CHKERRXX(ierr);
    } else {
      ierr = PetscViewerASCIIGetStdout(comm, &viewer); CHKERRXX(ierr);
    }
#if PETSC_VERSION_GE(3,7,0)
    ierr = PetscViewerPushFormat(viewer, format);
#else
    ierr = PetscViewerSetFormat(viewer, format);
#endif
    CHKERRXX(ierr);
    ierr = VecView(vec, viewer); CHKERRXX(ierr);
#if PETSC_VERSION_GE(3,7,0)
    ierr = PetscViewerPopFormat(viewer); CHKERRXX(ierr);
#endif
    if (filename != NULL) {
      ierr = PetscViewerDestroy(&viewer); CHKERRXX(ierr);
    }
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscVectorWrapper::print
// -------------------------------------------------------------
void
PetscVectorWrapper::print(const char* filename) const
{
  petsc_print_vector(p_vector, filename, PETSC_VIEWER_ASCII_INDEX);
}

// -------------------------------------------------------------
// PetscVectorWrapper::save
// -------------------------------------------------------------
void
PetscVectorWrapper::save(const char* filename) const
{
  petsc_print_vector(p_vector, filename, PETSC_VIEWER_ASCII_MATLAB);
}

// -------------------------------------------------------------
// PetscVectorWrapper::loadBinary
// -------------------------------------------------------------
void
PetscVectorWrapper::loadBinary(const char* filename)
{
  PetscErrorCode ierr;
  try {
    PetscViewer viewer;
    MPI_Comm comm = PetscObjectComm((PetscObject)p_vector);
    ierr = PetscViewerBinaryOpen(comm,
                                 filename,
                                 FILE_MODE_READ,
                                 &viewer); CHKERRXX(ierr);
    ierr = VecLoad(p_vector, viewer); CHKERRXX(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscVectorWrapper::saveBinary
// -------------------------------------------------------------
void
PetscVectorWrapper::saveBinary(const char* filename) const
{
  PetscErrorCode ierr;
  try {
    PetscViewer viewer;
    MPI_Comm comm = PetscObjectComm((PetscObject)p_vector);
    ierr = PetscViewerBinaryOpen(comm,
                                 filename,
                                 FILE_MODE_WRITE,
                                 &viewer); CHKERRXX(ierr);
    ierr = VecView(p_vector, viewer); CHKERRXX(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscVectorWrapper::p_accept
// -------------------------------------------------------------
void 
PetscVectorWrapper::p_accept(ImplementationVisitor& visitor)
{
  visitor.visit(*this);
}

void 
PetscVectorWrapper::p_accept(ConstImplementationVisitor& visitor) const
{
  visitor.visit(*this);
}



} // namespace math
} // namespace gridpack

