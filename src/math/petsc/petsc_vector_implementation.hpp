// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_vector_implementation.hpp
 * @author William A. Perkins
 * @date   2014-10-30 11:55:25 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_vector_implementation_h_
#define _petsc_vector_implementation_h_

#include "vector_implementation.hpp"
#include "petsc/petsc_vector_wrapper.hpp"
#include "petsc/petsc_exception.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScVectorImplementation
// -------------------------------------------------------------
/// Vector implementation based on the PETSc library
/**
 * 
 * 
 */
class PETScVectorImplementation 
  : public VectorImplementation {
public:

  /// Default constructor.
  /** 
   * @e Collective on @c comm.
   * 
   * @param comm parallel environment
   * @param local_length number of entries to be locally owned
   * 
   * @return 
   */
  PETScVectorImplementation(const parallel::Communicator& comm,
                            const IdxType& local_length)
    : VectorImplementation(comm), p_vwrap(comm, local_length)
  { }

  /// Construct from an existing PETSc vector
  PETScVectorImplementation(Vec& pvec, const bool& copyvec = true)
    : VectorImplementation(PetscVectorWrapper::getCommunicator(pvec)), 
      p_vwrap(pvec, copyvec)
  { }

  /// Destructor
  /** 
   * @e Collective
   * 
   */
  ~PETScVectorImplementation(void) {}

  // /// Get (a pointer to) the PETSc implementation

  // const Vec *getVector(void) const
  // {
  //   return p_vwrap.getVector();
  // }

  // /// Get (a pointer to) the PETSc implementation
  // Vec *getVector(void)
  // {
  //   return p_vwrap.getVector();
  // }

protected:
  
  /// Where the actual vector is stored
  PetscVectorWrapper p_vwrap;

  /// Get the global vector length
  IdxType p_size(void) const
  {
    return p_vwrap.size();
  }

  /// Get the size of the vector local part
  IdxType p_localSize(void) const
  {
    return p_vwrap.localSize();
  }

  /// Get the local min/max global indexes (specialized)
  void p_localIndexRange(IdxType& lo, IdxType& hi) const
  {
    PetscInt plo, phi;
    p_vwrap.localIndexRange(plo, phi);
    lo = plo;
    hi = phi;
  }

  /// Set an individual element (specialized)
  /** 
   * If you try to set an off processor value, it will be ignored
   * 
   * @param i global vector index
   * @param x value
   */
  void p_setElement(const IdxType& i, const TheType& x)
  {
    PetscErrorCode ierr;
    try {
      Vec *v = p_vwrap.getVector();
      ierr = VecSetValue(*v, i, x, INSERT_VALUES); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Set an several elements (specialized)
  void p_setElements(const IdxType& n, const IdxType *i, const TheType *x)
  {
    PetscErrorCode ierr;
    try {
      Vec *v = p_vwrap.getVector();
      ierr = VecSetValues(*v, n, i, x, INSERT_VALUES); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Add to an individual element (specialized)
  void p_addElement(const IdxType& i, const TheType& x)
  {
    PetscErrorCode ierr;
    try {
      Vec *v = p_vwrap.getVector();
      ierr = VecSetValue(*v, i, x, ADD_VALUES); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Add to an several elements (specialized)
  void p_addElements(const IdxType& n, const IdxType *i, const TheType *x)
  {
    PetscErrorCode ierr;
    try {
      Vec *v = p_vwrap.getVector();
      ierr = VecSetValues(*v, n, i, x, ADD_VALUES); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Get an individual (local) element (specialized)
  void p_getElement(const IdxType& i, TheType& x) const
  {
    this->p_getElements(1, &i, &x);
  }

  /// Get an several (local) elements (specialized)
  void p_getElements(const IdxType& n, const IdxType *i, TheType *x) const
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

  /// Get all of vector elements (on all processes)
  void p_getAllElements(TheType *x) const
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

  /// Make all the elements zero (specialized)
  void p_zero(void)
  {
    p_vwrap.zero();
  }

  /// Make all the elements the specified value (specialized)
  void p_fill(const TheType& v)
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

  /// Compute the vector L1 norm (sum of absolute value) (specialized)
  double p_norm1(void) const
  {
    return p_vwrap.norm1();
  }

  /// Compute the vector L2 norm (root of sum of squares) (specialized)
  double p_norm2(void) const
  {
    return p_vwrap.norm2();
  }

  /// Compute the vector infinity (or maximum) norm (specialized)
  double p_normInfinity(void) const
  {
    return p_vwrap.normInfinity();
  }

  /// Replace all elements with its absolute value (specialized) 
  void p_abs(void)
  {
    p_vwrap.abs();
  }

  /// Replace all elements with their complex conjugate
  void p_conjugate(void)
  {
    p_vwrap.conjugate();
  }

  /// Replace all elements with its exponential (specialized)
  void p_exp(void)
  {
    p_vwrap.exp();
  }

  /// Make this instance ready to use
  void p_ready(void)
  {
    p_vwrap.ready();
  }

  /// Allow visits by implemetation visitor
  void p_accept(ImplementationVisitor& visitor)
  {
    p_vwrap.accept(visitor);
  }

  /// Allow visits by implemetation visitor
  void p_accept(ConstImplementationVisitor& visitor) const
  {
    p_vwrap.accept(visitor);
  }

  /// Make an exact replica of this instance (specialized)
  VectorImplementation *p_clone(void) const
  {
    parallel::Communicator comm(this->communicator());
    IdxType local_size(this->localSize());
    
    PETScVectorImplementation *result = 
      new PETScVectorImplementation(comm, local_size);
    PetscErrorCode ierr;
    
    Vec *to_vec(result->p_vwrap.getVector());

    try {
      const Vec *v = p_vwrap.getVector();
      ierr = VecCopy(*v, *to_vec); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
    return result;
  }
};


} // namespace math
} // namespace gridpack

#endif
