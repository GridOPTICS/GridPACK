// -------------------------------------------------------------
/**
 * @file   petsc_vector_implementation.cpp
 * @author William A. Perkins
 * @date   2013-06-04 12:58:13 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#include "gridpack/math/implementation_visitor.hpp"
#include "gridpack/math/petsc/petsc_vector_implementation.hpp"
#include "gridpack/math/petsc/petsc_exception.hpp"
#include "gridpack/math/petsc/petsc_vector_extractor.hpp"

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
    ierr = VecCreate(comm,&p_vector); CHKERRXX(ierr);
    ierr = VecSetSizes(p_vector, local_length, PETSC_DECIDE); CHKERRXX(ierr);
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
// PETScVectorImplementation::p_local_size
// -------------------------------------------------------------
int 
PETScVectorImplementation::p_local_size(void) const
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
// PETScVectorImplementation::p_local_index_range
// -------------------------------------------------------------
void
PETScVectorImplementation::p_local_index_range(int& lo, int& hi) const
{
  lo = p_min_index;
  hi = p_max_index;
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_set_element
// -------------------------------------------------------------
/** 
 * If you try to set an off processor value, it will be ignored
 * 
 * @param i global vector index
 * @param x value
 */
void
PETScVectorImplementation::p_set_element(const int& i, const complex_type& x)
{
  PetscErrorCode ierr;
  try {
    ierr = VecSetValue(p_vector, i, x, INSERT_VALUES); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_set_elements
// -------------------------------------------------------------
void
PETScVectorImplementation::p_set_elements(const int& n, const int *i, const complex_type *x)
{
  PetscErrorCode ierr;
  try {
    ierr = VecSetValues(p_vector, n, i, x, INSERT_VALUES); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_add_element
// -------------------------------------------------------------
void
PETScVectorImplementation::p_add_element(const int& i, const complex_type& x)
{
  PetscErrorCode ierr;
  try {
    ierr = VecSetValue(p_vector, i, x, ADD_VALUES); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_add_elements
// -------------------------------------------------------------
void
PETScVectorImplementation::p_add_elements(const int& n, const int *i, const complex_type *x)
{
  PetscErrorCode ierr;
  try {
    ierr = VecSetValues(p_vector, n, i, x, ADD_VALUES); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_get_element
// -------------------------------------------------------------
void
PETScVectorImplementation::p_get_element(const int& i, complex_type& x) const
{
  this->p_get_elements(1, &i, &x);
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_get_elements
// -------------------------------------------------------------
void
PETScVectorImplementation::p_get_elements(const int& n, const int *i, complex_type *x) const
{
  PetscErrorCode ierr;
  try {
    ierr = VecGetValues(p_vector, n, i, x); CHKERRXX(ierr);
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
  complex_type v(0.0, 0.0);
  this->fill(v);
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_fill
// -------------------------------------------------------------
void
PETScVectorImplementation::p_fill(const complex_type& v)
{
  PetscErrorCode ierr(0);
  try {
    PetscScalar pv(v);
    ierr = VecSet(this->p_vector, pv);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
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
  int local_size(this->local_size());
  
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


// -------------------------------------------------------------
// PETScVectorImplementation::p_scale
// -------------------------------------------------------------
void 
PETScVectorImplementation::p_scale(const complex_type& x)
{
  PetscErrorCode ierr(0);
  try {
    ierr = VecScale(p_vector, x);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::add
// -------------------------------------------------------------
void
PETScVectorImplementation::p_add(const VectorImplementation& x)
{
  PetscErrorCode ierr(0);
  const Vec *xvec;
  {
    PETScConstVectorExtractor vext;
    x.accept(vext);
    xvec = vext.vector();
  }
  try {
    PetscScalar alpha(1.0);

    // This call computes y = x + alpha*y. Where y is p_vector.  
    ierr = VecAYPX(p_vector, alpha, *xvec);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

void
PETScVectorImplementation::p_add(const complex_type& x)
{
  PetscErrorCode ierr(0);
  try {
    ierr = VecShift(p_vector, x);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_copy
// -------------------------------------------------------------
void
PETScVectorImplementation::p_copy(const VectorImplementation& x)
{
  PetscErrorCode ierr(0);
  const Vec *xvec;
  {
    PETScConstVectorExtractor vext;
    x.accept(vext);
    xvec = vext.vector();
  }
  try {
    ierr = VecCopy(*xvec, p_vector);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScVectorImplementation::p_reciprocal
// -------------------------------------------------------------
void
PETScVectorImplementation::p_reciprocal(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = VecReciprocal(p_vector);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

} // namespace math
} // namespace gridpack

