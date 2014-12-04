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
 * @date   2014-10-28 13:41:59 d3g096
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
                            const IdxType& local_length);

  /// Construct from an existing PETSc vector
  PETScVectorImplementation(Vec& pvec, const bool& copyvec = true);

  /// Destructor
  /** 
   * @e Collective
   * 
   */
  ~PETScVectorImplementation(void);

  /// Get (a pointer to) the PETSc implementation

  const Vec *getVector(void) const;

  /// Get (a pointer to) the PETSc implementation
  Vec *getVector(void);

protected:
  
  /// Where the actual vector is stored
  PetscVectorWrapper p_vwrap;

  /// Get the global vector length
  IdxType p_size(void) const;

  /// Get the size of the vector local part
  IdxType p_localSize(void) const;

  /// Get the local min/max global indexes (specialized)
  void p_localIndexRange(IdxType& lo, IdxType& hi) const;

  /// Set an individual element (specialized)
  void p_setElement(const IdxType& i, const TheType& x);

  /// Set an several elements (specialized)
  void p_setElements(const IdxType& n, const IdxType *i, const TheType *x);

  /// Add to an individual element (specialized)
  void p_addElement(const IdxType& i, const TheType& x);

  /// Add to an several elements (specialized)
  void p_addElements(const IdxType& n, const IdxType *i, const TheType *x);

  /// Get an individual (local) element (specialized)
  void p_getElement(const IdxType& i, TheType& x) const;

  /// Get an several (local) elements (specialized)
  void p_getElements(const IdxType& n, const IdxType *i, TheType *x) const;

  /// Get all of vector elements (on all processes)
  void p_getAllElements(TheType *x) const;

  /// Make all the elements zero (specialized)
  void p_zero(void);

  /// Make all the elements the specified value (specialized)
  void p_fill(const TheType& v);

  /// Compute the vector L1 norm (sum of absolute value) (specialized)
  double p_norm1(void) const;

  /// Compute the vector L2 norm (root of sum of squares) (specialized)
  double p_norm2(void) const;

  /// Compute the vector infinity (or maximum) norm (specialized)
  double p_normInfinity(void) const;

  /// Replace all elements with its absolute value (specialized) 
  void p_abs(void);

  /// Replace all elements with their complex conjugate
  void p_conjugate(void);

  // FIXME: more ...

  /// Make this instance ready to use
  void p_ready(void);

  /// Allow visits by implemetation visitor
  void p_accept(ImplementationVisitor& visitor);

  /// Allow visits by implemetation visitor
  void p_accept(ConstImplementationVisitor& visitor) const;

  /// Make an exact replica of this instance (specialized)
  VectorImplementation *p_clone(void) const;

  // -------------------------------------------------------------
  // In-place Vector Operation Methods (change this instance)
  // -------------------------------------------------------------

};


} // namespace math
} // namespace gridpack

#endif
