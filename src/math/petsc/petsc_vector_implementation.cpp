// -------------------------------------------------------------
/**
 * @file   petsc_vector_implementation.cpp
 * @author William A. Perkins
 * @date   2013-05-08 14:06:13 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#include "gridpack/math/implementation_visitor.hpp"
#include "gridpack/math/petsc/petsc_vector_implementation.hpp"
#include "gridpack/math/petsc/petsc_exception.hpp"

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
  : VectorImplementation(comm)
{
  PetscErrorCode ierr;
  try {
    PetscInt lo, hi;
    ierr = VecCreate(comm,&vector_); CHKERRXX(ierr);
    ierr = VecSetSizes(vector_, local_length, PETSC_DECIDE); CHKERRXX(ierr);
    if (comm.size() > 1) {
      ierr = VecSetType(vector_, VECMPI);  CHKERRXX(ierr);
    } else {
      ierr = VecSetType(vector_, VECSEQ);  CHKERRXX(ierr);
    }
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
      ierr = VecDestroy(&vector_);
    }
  } catch (...) {
    // just eat it
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::size_
// -------------------------------------------------------------
int 
PETScVectorImplementation::size_(void) const
{
  PetscErrorCode ierr;
  try {
    PetscInt gsize;
    ierr = VecGetSize(this->vector_, &gsize); CHKERRXX(ierr);
    return static_cast<int>(gsize);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::local_size_
// -------------------------------------------------------------
int 
PETScVectorImplementation::local_size_(void) const
{
  PetscErrorCode ierr;
  try {
    PetscInt gsize;
    ierr = VecGetLocalSize(this->vector_, &gsize); CHKERRXX(ierr);
    return static_cast<int>(gsize);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::ready_
// -------------------------------------------------------------
void
PETScVectorImplementation::ready_(void)
{
  PetscErrorCode ierr;
  try {
    ierr = VecAssemblyBegin(vector_); CHKERRXX(ierr);
    ierr = VecAssemblyEnd(vector_); 
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }

}


// -------------------------------------------------------------
// PETScVectorImplementation::accept_
// -------------------------------------------------------------
void 
PETScVectorImplementation::accept_(ImplementationVisitor& visitor)
{
  visitor.visit(*this);
}


} // namespace math
} // namespace gridpack

