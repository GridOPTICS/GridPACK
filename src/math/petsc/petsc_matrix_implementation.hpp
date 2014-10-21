// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_implementation.h
 * @author William A. Perkins
 * @date   2014-10-21 12:55:19 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_matrix_implementation_h_
#define _petsc_matrix_implementation_h_

#include <petscmat.h>
#include "matrix_implementation.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScMatrixImplementation
// -------------------------------------------------------------
/// Matrix implementation based on the PETSc library
/**
 * 
 * 
 */
class PETScMatrixImplementation
  : public MatrixImplementation
{
public:

  /// Default constructor
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const IdxType& local_rows, const IdxType& local_cols,
                            const bool& dense = false);

  /// Construct a sparse matrix with an estimate of (maximum) usage
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const IdxType& local_rows, const IdxType& local_cols,
                            const IdxType& max_nonzero_per_row);

  /// Construct a sparse matrix with number of nonzeros in each row
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const IdxType& local_rows, const IdxType& local_cols,
                            const IdxType *nonzeros_by_row);
  

  /// Make a new instance from an existing PETSc matrix
  PETScMatrixImplementation(Mat& m, const bool& copymat = true);

  /// Destructor
  ~PETScMatrixImplementation(void);

  /// Get the PETSc matrix
  Mat *get_matrix(void)
  {
    return &p_matrix;
  }

  /// Get the PETSc matrix (const version)
  const Mat *get_matrix(void) const
  {
    return &p_matrix;
  }

protected:

  /// Extract a Communicator from a PETSc vector
  static parallel::Communicator p_getCommunicator(const Mat& m);

  /// The PETSc matrix representation
  Mat p_matrix;

  /// Was @c p_matrix created or just wrapped
  bool p_matrixWrapped;

  /// Build the generic PETSc matrix instance
  void p_build_matrix(const parallel::Communicator& comm,
                      const IdxType& local_rows, const IdxType& cols);

  /// Set up a dense matrix
  void p_set_dense_matrix(void);

  /// Set up a sparse matrix that is not preallocated
  void p_set_sparse_matrix(void);

  /// Set up a sparse matrix and preallocate it using the maximum nonzeros per row
  void p_set_sparse_matrix(const IdxType& max_nz_per_row);

  /// Set up a sparse matrix and preallocate it using known nonzeros for each row
  void p_set_sparse_matrix(const IdxType *nz_by_row);

  /// Get the global index range of the locally owned rows (specialized)
  void p_localRowRange(IdxType& lo, IdxType& hi) const;

  /// Get the total number of rows in this matrix (specialized)
  IdxType p_rows(void) const;

  /// Get the number of local rows in this matirx (specialized)
  IdxType p_localRows(void) const;

  /// Get the number of columns in this matrix (specialized)
  IdxType p_cols(void) const;

  /// Get the number of local rows in this matirx (specialized)
  IdxType p_localCols(void) const;

  /// Set an individual element
  void p_setElement(const IdxType& i, const IdxType& j, const TheType& x);

  /// Set an several element
  void p_setElements(const IdxType& n, const IdxType *i, const IdxType *j, const TheType *x);

  // /// Set all elements in a row
  // void p_set_row(const IdxType& i, const IdxType *j, const TheType *x);

  // /// Set all elements in a region
  // void p_set_region(const IdxType& ni, const IdxType& nj, 
  //                          const IdxType *i, const IdxType *j, const TheType *x);

  /// Add to  an individual element
  void p_addElement(const IdxType& i, const IdxType& j, const TheType& x);

  /// Add to  an several element
  void p_addElements(const IdxType& n, const IdxType *i, const IdxType *j, const TheType *x);

  // /// Add to  all elements in a row
  // void p_add_row(const IdxType& i, const IdxType *j, const TheType *x);

  /// Get an individual element
  void p_getElement(const IdxType& i, const IdxType& j, TheType& x) const;

  /// Get an several element
  void p_getElements(const IdxType& n, const IdxType *i, const IdxType *j, TheType *x) const;

  /// Replace all elements with their real parts
  void p_real(void);

  /// Replace all elements with their imaginary parts
  void p_imaginary(void);

  /// Replace all elements with their complex gradient
  void p_conjugate(void);

  /// Compute the matrix L<sup>2</sup> norm (specialized)
  double p_norm2(void) const;

  /// Make this instance ready to use
  void p_ready(void);

  /// Allow visits by implementation visitors
  void p_accept(ImplementationVisitor& visitor);

  /// Allow visits by implementation visitors
  void p_accept(ConstImplementationVisitor& visitor) const;

  /// Make an exact replica of this instance (specialized)
  MatrixImplementation *p_clone(void) const;

};

} // namespace utility
} // namespace gridpack



#endif
