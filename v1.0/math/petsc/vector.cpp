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
 * @date   2013-10-09 13:26:10 d3g096
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
  PETScVectorImplementation *impl = 
    new PETScVectorImplementation(comm, local_length);
  p_vector_impl.reset(impl);
  p_set_distributed(impl);
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
Vector::add(const Vector& x, const ComplexType& scale)
{
  this->p_checkCompatible(x);
  PetscErrorCode ierr(0);
  const Vec *xvec(PETScVector(x));
  Vec *yvec(PETScVector(*this));
  try {
    PetscScalar alpha(scale);

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
  this->p_checkCompatible(x);
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

// -------------------------------------------------------------
// petsc_make_viewer
// -------------------------------------------------------------
static void
petsc_make_viewer(const char* filename, PetscViewer *viewer)
{
  PetscErrorCode ierr;
  if (filename != NULL) {
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, viewer); ; CHKERRXX(ierr);
  } else {
    *viewer = PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD);
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
    petsc_make_viewer(filename, &viewer);
    ierr = PetscViewerSetFormat(viewer, format); ; CHKERRXX(ierr);
    ierr = VecView(vec, viewer); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}


// -------------------------------------------------------------
// Vector::print
// -------------------------------------------------------------
void
Vector::print(const char* filename) const
{
  const Vec *vec(PETScVector(*this));
  petsc_print_vector(*vec, filename, PETSC_VIEWER_ASCII_INDEX);
}

// -------------------------------------------------------------
// Vector::save
// -------------------------------------------------------------
void
Vector::save(const char* filename) const
{
  const Vec *vec(PETScVector(*this));
  petsc_print_vector(*vec, filename, PETSC_VIEWER_ASCII_MATLAB);
}


} // namespace math
} // namespace gridpack
