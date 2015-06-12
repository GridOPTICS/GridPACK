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
 * @date   2015-06-12 10:00:46 d3g096
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
#include "fallback_matrix_methods.hpp"


extern PetscErrorCode 
sillyMatScaleComplex(Mat A, const gridpack::ComplexType& px);

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
      p_mwrap()
  {
    IdxType tmp(max_nonzero_per_row*elementSize);
    p_mwrap.reset(new PetscMatrixWrapper(comm, 
                                         local_rows*elementSize, 
                                         local_cols*elementSize, 
                                         tmp));
  }

  /// Construct a sparse matrix with number of nonzeros in each row
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const IdxType& local_rows, const IdxType& local_cols,
                            const IdxType *nonzeros_by_row)
    : MatrixImplementation<T, I>(comm)
  {
    std::vector<IdxType> tmp(local_rows*elementSize);
    for (unsigned int i = 0; i < local_rows; ++i) {
      tmp[i*elementSize] = nonzeros_by_row[i]*elementSize;
      if (elementSize > 1) {
        tmp[i*elementSize+1] = nonzeros_by_row[i]*elementSize;
      }
    }
    p_mwrap.reset(new PetscMatrixWrapper(comm, 
                                         local_rows*elementSize, 
                                         local_cols*elementSize, 
                                         &tmp[0]));
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

  /// Create a sparse matrix with more or less specific ownership
  static 
  PETScMatrixImplementation *
  createDense(const parallel::Communicator& comm,
              const IdxType& global_rows, 
              const IdxType& global_cols,
              const IdxType& local_rows, 
              const IdxType& local_cols)
  {
    PetscErrorCode ierr(0);
    Mat mtmp;
    PETScMatrixImplementation *result;
    try {
      PetscInt grow(global_rows > 0 ? global_rows*elementSize : PETSC_DETERMINE);
      PetscInt gcol(global_cols > 0 ? global_cols*elementSize : PETSC_DETERMINE);
      PetscInt lrow(local_rows > 0 ? local_rows*elementSize : PETSC_DECIDE);
      PetscInt lcol(local_cols > 0 ? local_cols*elementSize : PETSC_DECIDE);

      ierr = MatCreate(comm, &mtmp); CHKERRXX(ierr);
      ierr = MatSetSizes(mtmp, lrow, lcol, grow, gcol); CHKERRXX(ierr);
      if (comm.size() == 1) {
        ierr = MatSetType(mtmp, MATSEQDENSE); CHKERRXX(ierr);
        ierr = MatSeqDenseSetPreallocation(mtmp, PETSC_NULL); CHKERRXX(ierr);
      } else {
        ierr = MatSetType(mtmp, MATDENSE); CHKERRXX(ierr);
        ierr = MatMPIDenseSetPreallocation(mtmp, PETSC_NULL); CHKERRXX(ierr);
      }
      result = new PETScMatrixImplementation(mtmp, false);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
    return result;
  }

  /// Get the PETSc matrix
  Mat *getMatrix(void)
  {
    return p_mwrap->getMatrix();
  }

  /// Get the PETSc matrix (const version)
  const Mat *getMatrix(void) const
  {
    return p_mwrap->getMatrix();
  }

protected:

  /// The actual PETSc matrix to be used
  boost::scoped_ptr<PetscMatrixWrapper> p_mwrap;

  /// Apply a specific unary operation to the vector
  void p_applyOperation(base_unary_function<TheType>& op)
  {
    PetscErrorCode ierr;

    Mat *A = this->getMatrix();
    PetscInt lo, hi;
    ierr = MatGetOwnershipRange(*A, &lo, &hi); CHKERRXX(ierr);

    boost::scoped_ptr< PETScMatrixImplementation<T, I> > 
      tmp(new PETScMatrixImplementation<T, I>(*A, true));

    std::vector<TheType> rvals;
    for (PetscInt i = lo; i < hi; i += 2) {
      PetscInt ncols;
      const PetscInt *cidx;
      const PetscScalar *p;
      int size;

      ierr = MatGetRow(*A, i, &ncols, &cidx, &p); CHKERRXX(ierr);

      // The arrays from MatGetRow are read only; make a copy of the values

      size = ncols/elementSize;
      rvals.resize(size);
      ValueTransferFromLibrary<PetscScalar, TheType> 
        trans(ncols, const_cast<PetscScalar *>(p), &rvals[0]);
      trans.go();
      
      // if TheType is complex and PetscScalar is real, the imaginary
      // part of all values will be negative; so copy the values to
      // another array and take the congjugate

      conjugate_value<TheType> c;
      std::transform(rvals.begin(), rvals.end(), rvals.begin(), c);

      // apply the operator

      for (typename std::vector<TheType>::iterator r = rvals.begin();
           r != rvals.end(); ++r) {
        *r = op(*r);
      }

      // put the new values in the temporary matrix
      for (int k = 0; k < ncols; k += elementSize) {
        int iidx(i/elementSize);
        int jidx(cidx[k]/elementSize);
        tmp->setElement(iidx, jidx, rvals[k/elementSize]);
      }

      ierr = MatRestoreRow(*A, i, &ncols, &cidx, &p); CHKERRXX(ierr);
    }
    tmp->ready();

    // copy the temporary matrix values back to the original matrix

    Mat *B = tmp->getMatrix();
    ierr = MatCopy(*B, *A, SAME_NONZERO_PATTERN); CHKERRXX(ierr);
  }

  /// Apply a specificy accumulator operation to the vector
  RealType p_applyAccumulator(base_accumulator_function<TheType, RealType>& op) const
  {
    RealType result(0.0);
    PetscErrorCode ierr;
    const Mat *A = this->getMatrix();
    
    PetscInt lo, hi;
    ierr = MatGetOwnershipRange(*A, &lo, &hi); CHKERRXX(ierr);
   
    std::vector<PetscScalar> rvals;
    for (PetscInt i = lo; i < hi; i += elementSize) {
      PetscInt ncols;
      const PetscInt *cidx;
      const PetscScalar *p;
      int size;

      ierr = MatGetRow(*A, i, &ncols, &cidx, &p); CHKERRXX(ierr);

      // The arrays from MatGetRow are read only; make a copy of the values

      size = ncols;
      rvals.resize(size);
      std::copy(p, p+size, rvals.begin());

      ierr = MatRestoreRow(*A, i, &ncols, &cidx, &p); CHKERRXX(ierr);

      // if TheType is complex and PetscScalar is real, the imaginary
      // part of all values will be negative; so copy the values to
      // another array and take the congjugate

      if (elementSize > 1) {
        conjugate_value<TheType> c;
        unary_operation<TheType, PetscScalar>(size, &rvals[0], c);
      }

      accumulator_operation<TheType, PetscScalar>(size, &rvals[0], op);

    }
    result = op.result();
    return result;
  }

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
    for (IdxType k = 0; k < n; k++) {
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
    for (IdxType k = 0; k < n; k++) {
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

          

  /// Scale this entire MatrixT by the given value (specialized)
  void p_scale(const TheType& xin)
  {
    PetscErrorCode ierr(0);
    Mat *pA(p_mwrap->getMatrix());
    if (useLibrary) {
      try {
        PetscScalar x = 
          gridpack::math::equate<PetscScalar, TheType>(xin);
        ierr = MatScale(*pA, x); CHKERRXX(ierr);
      } catch (const PETSC_EXCEPTION_TYPE& e) {
        throw PETScException(ierr, e);
      }
    } else {
      gridpack::math::multiplyvalue<TheType> op(xin);
      p_applyOperation(op);
    }
  }

  /// Shift the diagonal of this matrix by the specified value (specialized)
  void p_addDiagonal(const TheType& x)
  {
    if (useLibrary) {
      PetscScalar a = 
        gridpack::math::equate<PetscScalar, TheType>(x);
      Mat *pA(p_mwrap->getMatrix());
      PetscErrorCode ierr(0);
      try {
        ierr = MatShift(*pA, a); CHKERRXX(ierr);
      } catch (const PETSC_EXCEPTION_TYPE& e) {
        throw PETScException(ierr, e);
      }
    } else {
      fallback::addDiagonal(*this, x);
    }
  }

  /// Make this matrix the identity matrix (specialized)
  void p_identity(void)
  {
    Mat *pA(p_mwrap->getMatrix());

    PetscErrorCode ierr(0);
    try {
      PetscBool flag;
      PetscScalar one(1.0);
      ierr = MatAssembled(*pA, &flag); CHKERRXX(ierr);
      if (!flag) {
        int lo, hi;
        this->localRowRange(lo, hi);
        for (int i = lo; i < hi; ++i) {
          this->setElement(i, i, 1.0);
        }
        this->ready();
      } else {
        ierr = MatZeroEntries(*pA); CHKERRXX(ierr);
        ierr = MatShift(*pA, one); CHKERRXX(ierr);
      }
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  // ///  Get a row and put it in a local array (specialized)
  // void p_getRow(const IdxType& row, TheType *x) const
  // {
  //   PetscErrorCode ierr(0);
  //   try {
  //     Mat *mat = p_mwrap->getMatrix();
  //     MPI_Comm comm(PetscObjectComm((PetscObject)*mat));
  //     PetscInt n(this->cols()*elementSize);
      
  //     PetscInt ridx(row*elementSize);
  //     std::vector<PetscInt> cidx(n);

  //     IS irow, icol;
  //     ierr = ISCreateGeneral(comm, 1, &ridx, PETSC_COPY_VALUES, &irow); CHKERRXX(ierr);
  //     for (PetscInt j = 0; j < n; ++j) cidx[j] = j;
  //     ierr = ISCreateGeneral(comm, n, &cidx[0], PETSC_COPY_VALUES, &icol); CHKERRXX(ierr);

  //     Mat *sub;
  //     ierr = MatGetSubMatrices(*mat, 1, &irow, &icol, MAT_INITIAL_MATRIX, &sub); CHKERRXX(ierr);

  //     std::vector<PetscScalar> cval(n);
  //     ridx = 0;
  //     ierr = MatGetValues(*sub, 1, &ridx, n, &cidx[0], &cval[0]); CHKERRXX(ierr);
      
  //     MatrixValueTransferFromLibrary<PetscScalar, TheType> trans(n, &cval[0], x);
  //     trans.go();

  //     ierr = MatDestroyMatrices(1, &sub); CHKERRXX(ierr);

  //   } catch (const PETSC_EXCEPTION_TYPE& e) {
  //     throw PETScException(ierr, e);
  //   }
  // }

  /// Get some rows and put them in a local array (specialized)
  void p_getRowBlock(const IdxType& nrow, const IdxType *rows, TheType *x) const
  {
    PetscErrorCode ierr(0);
    try {
      Mat *mat = p_mwrap->getMatrix();
      MPI_Comm comm(PetscObjectComm((PetscObject)*mat));
      PetscInt ncol(this->cols()*elementSize);
      
      std::vector<PetscInt> ridx(nrow);
      std::vector<PetscInt> cidx(ncol);

      IS irow, icol;
      for (PetscInt i = 0; i < nrow; ++i) {
        ridx[i] = rows[i]*elementSize;
      }
      ierr = ISCreateGeneral(comm, nrow, &ridx[0], PETSC_COPY_VALUES, &irow); CHKERRXX(ierr);
      for (PetscInt j = 0; j < ncol; ++j) {
        cidx[j] = j;
      }
      ierr = ISCreateGeneral(comm, ncol, &cidx[0], PETSC_COPY_VALUES, &icol); CHKERRXX(ierr);

      Mat *sub;
      ierr = MatGetSubMatrices(*mat, 1, &irow, &icol, MAT_INITIAL_MATRIX, &sub); CHKERRXX(ierr);

      std::vector<PetscScalar> cval(nrow*ncol);
      for (PetscInt i = 0; i < nrow; ++i) {
        ridx[i] = i;
      }
      ierr = MatGetValues(*sub, nrow, &ridx[0], ncol, &cidx[0], &cval[0]); CHKERRXX(ierr);

      
      ValueTransferFromLibrary<PetscScalar, TheType> trans(nrow*ncol, &cval[0], x);
      trans.go();

      ierr = MatDestroyMatrices(1, &sub); CHKERRXX(ierr);

    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Replace all elements with their real parts
  void p_real(void)
  {
    if (useLibrary) {
      p_mwrap->real();
    } else {
      real_value<TheType> op;
      p_applyOperation(op);
    }
      
  }

  /// Replace all elements with their imaginary parts
  void p_imaginary(void)
  {
    if (useLibrary) {
      p_mwrap->imaginary();
    } else {
      imaginary_value<TheType> op;
      p_applyOperation(op);
    }
  }

  /// Replace all elements with their complex gradient
  void p_conjugate(void)
  {
    if (useLibrary) {
      p_mwrap->conjugate();
    } else {
      conjugate_value<TheType> op;
      p_applyOperation(op);
    }
  }

  /// Compute the matrix L<sup>2</sup> norm (specialized)
  double p_norm2(void) const
  {
    double result;
    if (useLibrary) {
      result = p_mwrap->norm2();
    } else {
      // l2_norm<TheType> op;
      // double lresult(p_applyAccumulator(op));
      // boost::mpi::all_reduce(this->communicator(), lresult, result, std::plus<double>());
      // result = sqrt(result);
      result = fallback::norm2<T, I>(this->communicator(), *this);
    }
    return result;
  }

  /// Zero all entries in the matrix (specialized)
  void p_zero(void)
  {
    p_mwrap->zero();
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
