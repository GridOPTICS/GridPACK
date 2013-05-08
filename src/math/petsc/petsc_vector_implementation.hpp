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
// Last Change: 2013-05-08 08:39:51 d3g096
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _petsc_vector_implementation_h_
#define _petsc_vector_implementation_h_

#include <petscvec.h>
#include "gridpack/math/vector_implementation.hpp"

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
    return &vector_;
  }

protected:

  /// The PETSc representation
  Vec vector_;

  /// Get the global vector length
  int size_(void) const;

  /// Get the size of the vector local part
  int local_size_(void) const;

  // /// Set an individual element (specialized)
  // void set_element_(const int& i, const complex_type& x);

  // /// Set an several elements (specialized)
  // void set_elements_(cont int& n, const int *i, const complex_type *x);

  // /// Add to an individual element (specialized)
  // void add_element_(const int& i, const complex_type& x);

  // /// Add to an several elements (specialized)
  // void add_elements_(const int& n, const int *i, const complex_type *x);

  // /// Make all the elements zero (specialized)
  // void zero_(void);

  // FIXME: more ...

  /// Make this instance ready to use
  void ready_(void);

  /// Allow visits by implemetation visitor
  void accept_(ImplementationVisitor& visitor);
};


} // namespace math
} // namespace gridpack

#endif
