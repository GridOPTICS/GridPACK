// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: petsc_vector_wrapper.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 24, 2014 by William A. Perkins
// Last Change: 2014-10-30 14:14:53 d3g096
// -------------------------------------------------------------


#ifndef _petsc_vector_wrapper_hpp_
#define _petsc_vector_wrapper_hpp_

#include <petscvec.h>
#include "parallel/communicator.hpp"
#include "implementation_visitable.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PetscVectorWrapper
// -------------------------------------------------------------
/**
 * Provide numeric type independent access to a PETSc Vec instance.
 * 
 */
class PetscVectorWrapper 
  : public ImplementationVisitable
{
public:

  /// Default constructor.
  PetscVectorWrapper(const parallel::Communicator& comm,
                     const PetscInt& local_length);

  /// Construct with an existing Vec instance
  PetscVectorWrapper(Vec& pVec, const bool& copyVec);

  /// Copy constructor
  PetscVectorWrapper(const PetscVectorWrapper& old);

  /// Destructor
  ~PetscVectorWrapper(void);

  /// Build a parallel::Communicator from the Vec instance
  static parallel::Communicator getCommunicator(const Vec& v);

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

  /// Get the global length
  /** 
   * 
   * @return global vector length
   */
  PetscInt size(void) const;

  /// Get the local length
  /** 
   * 
   * @return local vector length
   */
  PetscInt localSize(void) const;

  /// Get the local min/max global indexes
  /** 
   * @e Local.
   * 
   * The minimum index in the first global (0-based) index owned by
   * the local process. The maximum is one more than the last global
   * index owned. An example of usage:
   * 
     \code{.cpp}
     int lo, hi;
     my_vector.local_index_range(lo, hi);
     for (int i = lo; i < hi; ++i) {
       ComplexType x;
       x = ...;
       v.setElement(i, x);
     }
     \endcode
   * 
   * @param lo first (0-based) index of locally owned elements
   * @param hi one more than the last (0-based) index of locally owned elements
   */
  void localIndexRange(PetscInt& lo, PetscInt& hi) const;

  /// Make all the elements zero
  void zero(void);

  /// Compute the vector L1 norm (sum of absolute value)
  double norm1(void) const;

  /// Compute the vector L2 norm (root of sum of squares)
  double norm2(void) const;

  /// Compute the infinity (or maximum) norm
  double normInfinity(void) const;

  /// Replace all elements with its absolute value (complex magnitude) 
  void abs(void);

  /// Replace all elements with their complex conjugate
  void conjugate(void);

  /// Replace all elements with its exponential
  void exp(void);

  /// Replace all elements with its reciprocal
  void reciprocal(void);

  /// Make this instance ready to use
  void ready(void);

  /// Print to named file or standard output
  void print(const char* filename = NULL) const;

  /// Save, in MatLAB format, to named file (collective)
  void save(const char *filename) const;

  /// Load from a named file of whatever binary format the math library uses
  void loadBinary(const char *filename);

  /// Save to named file in whatever binary format the math library uses
  void saveBinary(const char *filename) const;

protected:

  /// Minimum global index on this processor
  PetscInt p_minIndex;

  /// Maximum global index on this processor
  PetscInt p_maxIndex;

  /// The PETSc representation
  Vec p_vector;

  /// Was @c p_vector created or just wrapped
  bool p_vectorWrapped;

  /// Compute Vec norm in the specified way
  double p_norm(const NormType& t) const;

  /// Allow visits by implemetation visitor
  void p_accept(ImplementationVisitor& visitor);

  /// Allow visits by implemetation visitor
  void p_accept(ConstImplementationVisitor& visitor) const;

};

} // namespace math
} // namespace gridpack

#endif
