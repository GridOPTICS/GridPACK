// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: petsc_vector_implementation.h
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 26, 2013 by William A. Perkins
// Last Change: 2013-06-11 14:13:06 d3g096
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _petsc_vector_implementation_h_
#define _petsc_vector_implementation_h_

#include <petscvec.h>
#include "vector_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScVectorImplementation
// -------------------------------------------------------------
class PETScVectorImplementation 
  : public VectorImplementation {
public:

  /// Default constructor.
  PETScVectorImplementation(const parallel::Communicator& comm,
                            const int& local_length);

  /// Destructor
  ~PETScVectorImplementation(void);

  /// Get (a pointer to) the PETSc implementation
  const Vec *get_vector(void) const
  {
    return &p_vector;
  }

  /// Get (a pointer to) the PETSc implementation
  Vec *get_vector(void)
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
  int p_local_size(void) const;

  /// Get the local min/max global indexes (specialized)
  void p_local_index_range(int& lo, int& hi) const;

  /// Set an individual element (specialized)
  void p_set_element(const int& i, const ComplexType& x);

  /// Set an several elements (specialized)
  void p_set_elements(const int& n, const int *i, const ComplexType *x);

  /// Add to an individual element (specialized)
  void p_add_element(const int& i, const ComplexType& x);

  /// Add to an several elements (specialized)
  void p_add_elements(const int& n, const int *i, const ComplexType *x);

  /// Get an individual element (specialized)
  void p_get_element(const int& i, ComplexType& x) const;

  /// Get an several elements (specialized)
  void p_get_elements(const int& n, const int *i, ComplexType *x) const;

  /// Make all the elements zero (specialized)
  void p_zero(void);

  /// Make all the elements the specified value (specialized)
  void p_fill(const ComplexType& v);

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
