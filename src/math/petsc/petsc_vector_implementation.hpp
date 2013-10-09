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
 * @date   2013-10-09 13:25:58 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_vector_implementation_h_
#define _petsc_vector_implementation_h_

#include <petscvec.h>
#include "vector_implementation.hpp"

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
                            const int& local_length);

  /// Destructor
  /** 
   * @e Collective
   * 
   */
  ~PETScVectorImplementation(void);

  /// Get (a pointer to) the PETSc implementation

  const Vec *getVector(void) const
  {
    return &p_vector;
  }

  /// Get (a pointer to) the PETSc implementation
  Vec *getVector(void)
  {
    return &p_vector;
  }

protected:

  /// Minimum global index on this processor
  int p_min_index;

  /// Maximum global index on this processor
  int p_max_index;

  /// The PETSc representation
  Vec p_vector;

  /// Get the global vector length
  int p_size(void) const;

  /// Get the size of the vector local part
  int p_localSize(void) const;

  /// Get the local min/max global indexes (specialized)
  void p_localIndexRange(int& lo, int& hi) const;

  /// Set an individual element (specialized)
  void p_setElement(const int& i, const ComplexType& x);

  /// Set an several elements (specialized)
  void p_setElements(const int& n, const int *i, const ComplexType *x);

  /// Add to an individual element (specialized)
  void p_addElement(const int& i, const ComplexType& x);

  /// Add to an several elements (specialized)
  void p_addElements(const int& n, const int *i, const ComplexType *x);

  /// Get an individual (local) element (specialized)
  void p_getElement(const int& i, ComplexType& x) const;

  /// Get an several (local) elements (specialized)
  void p_getElements(const int& n, const int *i, ComplexType *x) const;

  /// Get all of vector elements (on all processes)
  void p_getAllElements(ComplexType *x) const;

  /// Make all the elements zero (specialized)
  void p_zero(void);

  /// Make all the elements the specified value (specialized)
  void p_fill(const ComplexType& v);

  /// Compute the vector L1 norm (sum of absolute value) (specialized)
  ComplexType p_norm1(void) const;

  /// Compute the vector L2 norm (root of sum of squares) (specialized)
  ComplexType p_norm2(void) const;

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
