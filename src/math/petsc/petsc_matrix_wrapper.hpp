// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: petsc_matrix_wrapper.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created January 28, 2015 by William A. Perkins
// Last Change: 2019-05-08 13:53:18 d3g096
// -------------------------------------------------------------


#ifndef _petsc_matrix_wrapper_hpp_
#define _petsc_matrix_wrapper_hpp_

#include <petscmat.h>
#include "parallel/communicator.hpp"
#include "implementation_visitable.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PetscMatrixWrapper
// -------------------------------------------------------------
/**
 * A numeric type independent wrapper for PETSc Mat instances
 * 
 */
class PetscMatrixWrapper 
  : public ImplementationVisitable
{
public:

  /// Extract a Communicator from a PETSc vector
  static parallel::Communicator getCommunicator(const Mat& m);

  /// Default constructor.
  PetscMatrixWrapper(const parallel::Communicator& comm,
                     const PetscInt& local_rows, const PetscInt& local_cols,
                     const bool& dense = false);

  /// Construct a sparse matrix allocating the same number of nonzeros in all rows
  PetscMatrixWrapper(const parallel::Communicator& comm,
                     const PetscInt& local_rows, const PetscInt& local_cols,
                     const PetscInt& max_nonzero_per_row);

  /// Construct a sparse matrix with nonzero count specified for each (local) row
  PetscMatrixWrapper(const parallel::Communicator& comm,
                     const PetscInt& local_rows, const PetscInt& local_cols,
                     const PetscInt *nonzeros_by_row);

  /// Constructor that wraps an existing Mat instance
  PetscMatrixWrapper(Mat& m, const bool& copymat = true,
                     const bool& destroymat = false);

  /// Destructor
  ~PetscMatrixWrapper(void);

  /// Get the total number of rows in this matrix (specialized)
  PetscInt rows(void) const;

  /// Get the number of local rows in this matirx (specialized)
  PetscInt localRows(void) const;

  /// Get the global index range of the locally owned rows (specialized)
  void localRowRange(PetscInt& lo, PetscInt& hi) const;

  /// Get the number of columns in this matrix (specialized)
  PetscInt cols(void) const;

  /// Get the number of local rows in this matirx (specialized)
  PetscInt localCols(void) const;

  /// Replace all elements with their real parts
  void real(void);

  /// Replace all elements with their imaginary parts
  void imaginary(void);

  /// Replace all elements with their complex gradient
  void conjugate(void);

  /// Compute the matrix L<sup>2</sup> norm (specialized)
  double norm2(void) const;

  /// Zero all entries in the matrix
  void zero(void);

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

  /// Get the PETSc matrix
  Mat *getMatrix(void)
  {
    return &p_matrix;
  }

  /// Get the PETSc matrix (const version)
  const Mat *getMatrix(void) const
  {
    return &p_matrix;
  }

protected:

  /// The real PETSc matrix
  Mat p_matrix;

  /// Was @c p_matrix created or just wrapped
  bool p_matrixWrapped;

  /// Destroy wrapped @c p_matrix even if it's wrapped
  bool p_destroyWrapped;

  /// Build the generic PETSc matrix instance
  void p_build_matrix(const parallel::Communicator& comm,
                      const PetscInt& local_rows, const PetscInt& cols);

  /// Set up a dense matrix
  void p_set_dense_matrix(void);

  /// Set up a sparse matrix that is not preallocated
  void p_set_sparse_matrix(void);

  /// Set up a sparse matrix and preallocate it using the maximum nonzeros per row
  void p_set_sparse_matrix(const PetscInt& max_nz_per_row);

  /// Set up a sparse matrix and preallocate it using known nonzeros for each row
  void p_set_sparse_matrix(const PetscInt *nz_by_row);

  /// Allow visits by implemetation visitor
  void p_accept(ImplementationVisitor& visitor);

  /// Allow visits by implemetation visitor
  void p_accept(ConstImplementationVisitor& visitor) const;

};

} // namespace math
} // namespace gridpack

#endif
