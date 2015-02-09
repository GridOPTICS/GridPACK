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
 * @date   2015-01-28 14:23:58 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_matrix_implementation_h_
#define _petsc_matrix_implementation_h_

#include <petscmat.h>
#include "petsc_exception.hpp"
#include "matrix_implementation.hpp"
#include "petsc_matrix_wrapper.hpp"

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
template <typename T, typename I = int>
class PETScMatrixImplementation
  : public MatrixImplementation<T, I>
{
public:

  typedef typename MatrixImplementation<T, I>::TheType TheType;
  typedef typename MatrixImplementation<T, I>::IdxType IdxType;

  /// Default constructor
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const IdxType& local_rows, const IdxType& local_cols,
                            const bool& dense = false)
    : MatrixImplementation<T, I>(comm),
      p_mwrap(comm, local_rows, local_cols, dense)
  {
  }

  /// Construct a sparse matrix with an estimate of (maximum) usage
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const IdxType& local_rows, const IdxType& local_cols,
                            const IdxType& max_nonzero_per_row)
    : MatrixImplementation<T, I>(comm),
      p_mwrap(comm, local_rows, local_cols, max_nonzero_per_row)
  {
  }

  /// Construct a sparse matrix with number of nonzeros in each row
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const IdxType& local_rows, const IdxType& local_cols,
                            const IdxType *nonzeros_by_row)
    : MatrixImplementation<T, I>(comm),
      p_mwrap(comm, local_rows, local_cols, nonzeros_by_row)
  {
  }

  /// Make a new instance from an existing PETSc matrix
  PETScMatrixImplementation(Mat& m, const bool& copyMat = true)
    : MatrixImplementation<T, I>(PetscMatrixWrapper::getCommunicator(m)),
      p_mwrap(m, copyMat)
  {
  }

  /// Destructor
  ~PETScMatrixImplementation(void)
  {
  }

  /// Get the PETSc matrix
  Mat *get_matrix(void)
  {
    return p_mwrap.getMatrix();
  }

  /// Get the PETSc matrix (const version)
  const Mat *get_matrix(void) const
  {
    return p_mwrap.getMatrix();
  }

protected:

  /// The actual PETSc matrix to be used
  PetscMatrixWrapper p_mwrap;

  /// Get the global index range of the locally owned rows (specialized)
  void p_localRowRange(IdxType& lo, IdxType& hi) const
  {
    PetscInt plo, phi;
    p_mwrap.localRowRange(plo, phi);
    lo = plo;
    hi = phi;
  }

  /// Get the total number of rows in this matrix (specialized)
  IdxType p_rows(void) const
  {
    return p_mwrap.rows();
  }

  /// Get the number of local rows in this matirx (specialized)
  IdxType p_localRows(void) const
  {
    return p_mwrap.localRows();
  }

  /// Get the number of columns in this matrix (specialized)
  IdxType p_cols(void) const
  {
    return p_mwrap.cols();
  }

  /// Get the number of local rows in this matirx (specialized)
  IdxType p_localCols(void) const
  {
    return p_mwrap.localCols();
  }

  /// Set an individual element
  void p_setElement(const IdxType& i, const IdxType& j, const TheType& x)
  {
    PetscErrorCode ierr(0);
    try {
      Mat *mat = p_mwrap.getMatrix();
      ierr = MatSetValue(*mat, i, j, x, INSERT_VALUES); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Set an several element
  void p_setElements(const IdxType& n, const IdxType *i, const IdxType *j, const TheType *x)
  {
    // FIXME: There's probably a better way
    for (PETScMatrixImplementation::IdxType k = 0; k < n; k++) {
      this->p_setElement(i[k], j[k], x[k]);
    }
  }

  /// Add to  an individual element
  void p_addElement(const IdxType& i, const IdxType& j, const TheType& x)
  {
    PetscErrorCode ierr(0);
    try {
      Mat *mat = p_mwrap.getMatrix();
      ierr = MatSetValue(*mat, i, j, x, ADD_VALUES); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Add to  an several element
  void p_addElements(const IdxType& n, const IdxType *i, const IdxType *j, const TheType *x)
  {
    // FIXME: There's probably a better way
    for (PETScMatrixImplementation::IdxType k = 0; k < n; k++) {
      this->p_addElement(i[k], j[k], x[k]);
    }
  }

  /// Get an individual element
  void p_getElement(const IdxType& i, const IdxType& j, TheType& x) const
  {
    PetscErrorCode ierr(0);
    try {
      const Mat *mat = p_mwrap.getMatrix();
      static const PETScMatrixImplementation::IdxType one(1);
      PetscInt iidx[one] = { i };
      PetscInt jidx[one] = { j };
      ierr = MatGetValues(*mat, one, &iidx[0], one, &jidx[0], &x); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Get an several element
  void p_getElements(const IdxType& n, const IdxType *i, const IdxType *j, TheType *x) const
  {
    // FIXME: There is a better way
    for (int k = 0; k < n; k++) {
      this->p_getElement(i[k], j[k], x[k]);
    }
  }

  /// Replace all elements with their real parts
  void p_real(void)
  {
    p_mwrap.real();
  }

  /// Replace all elements with their imaginary parts
  void p_imaginary(void)
  {
    p_mwrap.imaginary();
  }

  /// Replace all elements with their complex gradient
  void p_conjugate(void)
  {
    p_mwrap.conjugate();
  }

  /// Compute the matrix L<sup>2</sup> norm (specialized)
  double p_norm2(void) const
  {
    return p_mwrap.norm2();
  }

  /// Make this instance ready to use
  void p_ready(void)
  {
    p_mwrap.ready();
  }

  /// Print to named file or standard output
  void p_print(const char* filename = NULL) const
  {
    p_mwrap.print(filename);
  }

  /// Save, in MatLAB format, to named file (collective)
  void p_save(const char *filename) const
  {
    p_mwrap.save(filename);
  }

  /// Load from a named file of whatever binary format the math library uses
  void p_loadBinary(const char *filename)
  {
    p_mwrap.loadBinary(filename);
  }

  /// Save to named file in whatever binary format the math library uses
  void p_saveBinary(const char *filename) const
  {
    p_mwrap.saveBinary(filename);
  }

  /// Allow visits by implementation visitors
  void p_accept(ImplementationVisitor& visitor)
  {
    p_mwrap.accept(visitor);
  }

  /// Allow visits by implementation visitors
  void p_accept(ConstImplementationVisitor& visitor) const
  {
    p_mwrap.accept(visitor);
  }

  /// Make an exact replica of this instance (specialized)
  MatrixImplementation<T, I> *p_clone(void) const
  {
    PETScMatrixImplementation<T, I> *result =
      new PETScMatrixImplementation<T, I>(*(const_cast<Mat*>(p_mwrap.getMatrix())), true);
    return result;
  }


};

} // namespace utility
} // namespace gridpack



#endif
