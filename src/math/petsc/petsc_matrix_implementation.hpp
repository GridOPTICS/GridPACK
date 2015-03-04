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
 * @date   2015-02-25 14:23:54 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_matrix_implementation_h_
#define _petsc_matrix_implementation_h_

#include <petscmat.h>
#include <boost/scoped_ptr.hpp>
#include "petsc_exception.hpp"
#include "petsc_types.hpp"
#include "matrix_implementation.hpp"
#include "petsc_matrix_wrapper.hpp"
#include "value_transfer.hpp"

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

  /// A flag to denote whether the library can be used directly
  /**
   * Some operations can be passed directly to the underlying library
   * if the TheType is the same as the PETSc type @e or the PETSc type
   * is complex.  This type computes and stores that flag. 
   * 
   */
  static const bool useLibrary = UsePetscLibrary<TheType>::value;

  /// The number of library elements used to represent a single vector element
  static const unsigned int elementSize = PetscElementSize<TheType>::value;
  

  /// Default constructor
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const IdxType& local_rows, const IdxType& local_cols,
                            const bool& dense = false)
    : MatrixImplementation<T, I>(comm),
      p_mwrap(new PetscMatrixWrapper(comm, 
                                     local_rows*elementSize, 
                                     local_cols*elementSize, 
                                     dense))
  {
  }

  /// Construct a sparse matrix with an estimate of (maximum) usage
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const IdxType& local_rows, const IdxType& local_cols,
                            const IdxType& max_nonzero_per_row)
    : MatrixImplementation<T, I>(comm),
      p_mwrap(new PetscMatrixWrapper(comm, 
                                     local_rows*elementSize, 
                                     local_cols*elementSize, 
                                     max_nonzero_per_row))
  {
  }

  /// Construct a sparse matrix with number of nonzeros in each row
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const IdxType& local_rows, const IdxType& local_cols,
                            const IdxType *nonzeros_by_row)
    : MatrixImplementation<T, I>(comm)
  {
    if (elementSize > 1) {
      std::vector<IdxType> tmp(local_rows*elementSize);
      for (unsigned int i = 0; i < tmp.size(); i += elementSize) {
        tmp[i] = nonzeros_by_row[i];
        tmp[i+1] = nonzeros_by_row[i];
      }
      p_mwrap.reset(new PetscMatrixWrapper(comm, 
                                           local_rows*elementSize, 
                                           local_cols*elementSize, 
                                           &tmp[0]));
    } else {
      p_mwrap.reset(new PetscMatrixWrapper(comm, 
                                           local_rows*elementSize, 
                                           local_cols*elementSize, 
                                           nonzeros_by_row));
    }
                    
  }

  /// Make a new instance from an existing PETSc matrix
  PETScMatrixImplementation(Mat& m, const bool& copyMat = true)
    : MatrixImplementation<T, I>(PetscMatrixWrapper::getCommunicator(m)),
      p_mwrap(new PetscMatrixWrapper(m, copyMat))
  {
  }

  /// Destructor
  ~PETScMatrixImplementation(void)
  {
  }

  /// Get the PETSc matrix
  Mat *get_matrix(void)
  {
    return p_mwrap->getMatrix();
  }

  /// Get the PETSc matrix (const version)
  const Mat *get_matrix(void) const
  {
    return p_mwrap->getMatrix();
  }

protected:

  /// The actual PETSc matrix to be used
  boost::scoped_ptr<PetscMatrixWrapper> p_mwrap;

  /// Get the global index range of the locally owned rows (specialized)
  void p_localRowRange(IdxType& lo, IdxType& hi) const
  {
    PetscInt plo, phi;
    p_mwrap->localRowRange(plo, phi);
    lo = plo/elementSize;
    hi = phi/elementSize;
  }

  /// Get the total number of rows in this matrix (specialized)
  IdxType p_rows(void) const
  {
    return p_mwrap->rows()/elementSize;
  }

  /// Get the number of local rows in this matirx (specialized)
  IdxType p_localRows(void) const
  {
    return p_mwrap->localRows()/elementSize;
  }

  /// Get the number of columns in this matrix (specialized)
  IdxType p_cols(void) const
  {
    return p_mwrap->cols()/elementSize;
  }

  /// Get the number of local rows in this matirx (specialized)
  IdxType p_localCols(void) const
  {
    return p_mwrap->localCols()/elementSize;
  }

  /// Set an individual element
  void p_setElement(const IdxType& i, const IdxType& j, const TheType& x, 
                    InsertMode mode)
  {
    PetscErrorCode ierr(0);
    try {
      Mat *mat = p_mwrap->getMatrix();
      TheType tmp(x);
      PetscScalar px[elementSize*elementSize];
      MatrixValueTransferToLibrary<TheType, PetscScalar> trans(1, &tmp, &px[0]);
      trans.go();
      int n(elementSize);
      PetscInt iidx[elementSize], jidx[elementSize];
      for (int ii = 0; ii < elementSize; ++ii) {
        iidx[ii] = i*elementSize + ii;
      }
      for (int jj = 0; jj < elementSize; ++jj) {
        jidx[jj] = j*elementSize + jj;
      }
      ierr = MatSetValues(*mat, n, &iidx[0], n, &jidx[0], &px[0], mode); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }


  void p_setElement(const IdxType& i, const IdxType& j, const TheType& x)
  {
    p_setElement(i, j, x, INSERT_VALUES);
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
    p_setElement(i, j, x, ADD_VALUES);
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
      Mat *mat = p_mwrap->getMatrix();
      PetscScalar px[elementSize*elementSize];
      int n(elementSize);
      PetscInt iidx[elementSize], jidx[elementSize];
      for (int ii = 0; ii < elementSize; ++ii) {
        iidx[ii] = i*elementSize + ii;
      }
      for (int jj = 0; jj < elementSize; ++jj) {
        jidx[jj] = j*elementSize + jj;
      }
      ierr = MatGetValues(*mat, n, &iidx[0], n, &jidx[0], &px[0]); CHKERRXX(ierr);
      MatrixValueTransferFromLibrary<PetscScalar, TheType> trans(1, &px[0], &x);
      trans.go();
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
    p_mwrap->real();
  }

  /// Replace all elements with their imaginary parts
  void p_imaginary(void)
  {
    p_mwrap->imaginary();
  }

  /// Replace all elements with their complex gradient
  void p_conjugate(void)
  {
    p_mwrap->conjugate();
  }

  /// Compute the matrix L<sup>2</sup> norm (specialized)
  double p_norm2(void) const
  {
    return p_mwrap->norm2();
  }

  /// Make this instance ready to use
  void p_ready(void)
  {
    p_mwrap->ready();
  }

  /// Print to named file or standard output
  void p_print(const char* filename = NULL) const
  {
    p_mwrap->print(filename);
  }

  /// Save, in MatLAB format, to named file (collective)
  void p_save(const char *filename) const
  {
    p_mwrap->save(filename);
  }

  /// Load from a named file of whatever binary format the math library uses
  void p_loadBinary(const char *filename)
  {
    p_mwrap->loadBinary(filename);
  }

  /// Save to named file in whatever binary format the math library uses
  void p_saveBinary(const char *filename) const
  {
    p_mwrap->saveBinary(filename);
  }

  /// Allow visits by implementation visitors
  void p_accept(ImplementationVisitor& visitor)
  {
    p_mwrap->accept(visitor);
  }

  /// Allow visits by implementation visitors
  void p_accept(ConstImplementationVisitor& visitor) const
  {
    p_mwrap->accept(visitor);
  }

  /// Make an exact replica of this instance (specialized)
  MatrixImplementation<T, I> *p_clone(void) const
  {
    PETScMatrixImplementation<T, I> *result =
      new PETScMatrixImplementation<T, I>(*(const_cast<Mat*>(p_mwrap->getMatrix())), true);
    return result;
  }


};

} // namespace utility
} // namespace gridpack



#endif
