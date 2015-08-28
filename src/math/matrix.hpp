// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   matrix.hpp
 * @author William A. Perkins
 * @date   2015-08-18 14:09:32 d3g096
 * 
 * @brief  Declaration of the Matrix class.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 22, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _matrix_hpp_
#define _matrix_hpp_

#include <boost/scoped_ptr.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/math/matrix_implementation.hpp>
#include <gridpack/math/vector.hpp>
#include <gridpack/math/matrix_storage_type.hpp>


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class Matrix
// -------------------------------------------------------------
/// A parallel or serial matrix of real values
/**
 * This class uses the Pimpl idiom for implementation in order so the
 * interface is completely free of the underlying library.  If
 * constructed with a parallel environment with only one process, a
 * serial storage scheme is created, otherwise it's parallel. 
 *
 */
template <typename T, typename I = int>
class MatrixT
  : public parallel::WrappedDistributed,
    private utility::Uncopyable,
    public BaseMatrixInterface<T, I>
{
public:

  typedef typename BaseMatrixInterface<T, I>::TheType TheType;
  typedef typename BaseMatrixInterface<T, I>::IdxType IdxType;

  /// Constructor.
  /** 
   * A Matrix must be instantiated simulutaneously on all processes
   * involved in the specified \ref parallel::Communicator
   * "communicator". Each process in the communicator will own the
   * number of rows requested.
   * 
   * @param dist parallel environment
   * @param local_rows matrix rows to be owned by the local process
   * @param local_cols matrix columns to be owned by the local process
   * @param storage_type specify dense or sparse storage 
   * 
   * @return 
   */
  MatrixT(const parallel::Communicator& dist,
          const int& local_rows,
          const int& local_cols,
          const MatrixStorageType& storage_type = Sparse);

  /// Sparse matrix constructor with maximum number of nonzeros in a row
  /** 
   * If the underlying math implementation supports it, this
   * constructs a sparse matrix and pre-allocates it to allow @c
   * max_nz_per_row nonzeros in each row. If the underlying math
   * library supports it, @c max_nz_per_row does not need to be the
   * same on all processors.
   * 
   * @param dist parallel environment
   * @param local_rows matrix rows to be owned by the local process
   * @param local_cols matrix columns to be owned by the local process
   * @param max_nz_per_row maximum number of nonzeros in a row
   * 
   * @return new MatrixT
   */
  MatrixT(const parallel::Communicator& dist,
          const int& local_rows,
          const int& local_cols,
          const int& max_nz_per_row);
  
  /// Sparse matrix constructor with number of nonzeros for each row
  /** 
   * 
   * 
   * @param dist parallel environment
   * @param local_rows matrix rows to be owned by the local process
   * @param local_cols matrix columns to be owned by the local process
   * @param nz_by_row 
   * 
   * @return 
   */
  MatrixT(const parallel::Communicator& dist,
          const int& local_rows,
          const int& local_cols,
          const int *nz_by_row);

  /// Construct with an existing (allocated) implementation 
  /** 
   * For internal use only.
   * 
   * @param impl 
   * 
   * @return 
   */
  MatrixT(MatrixImplementation<T, I> *impl)
    : parallel::WrappedDistributed(impl), 
      utility::Uncopyable(),
      p_matrix_impl(impl)
  {
    BOOST_ASSERT(p_matrix_impl);
  }

  /// Destructor
  /** 
   * A matrix must be destroyed simulutaneously on all processes
   * involved in the \ref parallel::Communicator "communicator" used
   * to instantiate it.
   */
  ~MatrixT(void)
  {
  }

  /// Create a ::Dense Matrix instance with more/less specific ownership
  /** 
   * 
   * 
   * @param comm 
   * @param stype 
   * @param global_rows 
   * @param global_cols 
   * @param local_rows 
   * @param local_cols 
   * 
   * @return 
   */
  static MatrixT *
  createDense(const parallel::Communicator& comm,
              const int& global_rows,
              const int& global_cols,
              const int& local_rows,
              const int& local_cols);
  
  /// Get the storage type of this matrix
  MatrixStorageType storageType(void) const;

  //! @cond DEVDOC

  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor)
  {
    p_matrix_impl->accept(visitor);
  }

  /// Allow visits by implemetation visitor
  void accept(ConstImplementationVisitor& visitor) const
  {
    p_matrix_impl->accept(visitor);
  }

  //! @endcond

  /// Make an exact replica of this instance
  /** 
   * @e Collective.
   *
   * 
   * 
   * 
   * @return pointer to new (allocated) matrix instance
   */
  MatrixT *clone(void) const
  {
    MatrixImplementation<T, I> *pimpl_clone =
      this->p_matrix_impl->clone();
    MatrixT *result = new MatrixT(pimpl_clone);
    return result;
  }

  MatrixT *localClone(void) const
  {
    MatrixImplementation<T, I> *pimpl_clone =
      this->p_matrix_impl->localClone();
    MatrixT *result = new MatrixT(pimpl_clone);
    return result;
  }

  // -------------------------------------------------------------
  // Matrix Operation Methods
  // -------------------------------------------------------------

  /// Get a part of the matrix (new instance allocated)
  // Matrix *submatrix(const int& istart, const int& jstart,
  //                   const int& iend, const int& jend) const;

  // -------------------------------------------------------------
  // In-Place Matrix Operation Methods (change this instance)
  // -------------------------------------------------------------

  /// Make this Matrix equal to another
  /** 
   * @e Collective.
   * 
   * @param A 
   */
  void equate(const MatrixT& A);

  /// Multiply this matrix diagonal by the specified vector
  /** 
   * @e Collective.
   *
   * This is element by element multiplication
   * 
   * @param x factor by which all diagonal elements in the matrix are multiplied
   */
  void multiplyDiagonal(const VectorT<T, I>& x);

  /// Add the specified vector to the diagonal of this matrix
  /** 
   * @c Collective.
   * 
   * @param x 
   */
  void addDiagonalVector(const VectorT<T, I>& x);

  /// Add another matrix to this one, in place
  /** 
   * @e Collective.
   *
   * The specified matrix @c A must be the same global size (rows and
   * columns) as this instance, but local ownership is not important,
   * nor are any differences in nonzero entry patterns.
   * 
   * @param A matrix to add to this instance
   */
  void add(const MatrixT& A);

  // -------------------------------------------------------------
  // MatrixT Operations 
  //
  // all allocate new instances and throw when arguments are
  // inconsistent
  // -------------------------------------------------------------
  // friend MatrixT *factorize(const MatrixT& A);
  // friend MatrixT *inverse(const MatrixT& A);
  // friend MatrixT *reorder(const MatrixT& A, const Reordering& r);
  // friend MatrixT *identity(const MatrixT& A);
  // friend template <typename T, typename I>
  // MatrixT<T, I> *multiply(const MatrixT<T, I>& A, const MatrixT<T, I>& B);
  // friend template <typename T, typename I>
  // MatrixT<T, I> *transpose(const MatrixT<T, I>& A);

protected:

  /// The actual implementation
  boost::scoped_ptr< MatrixImplementation<T, I> > p_matrix_impl;

  /// Get the global index range of the locally owned rows (specialized)
  void p_localRowRange(IdxType& lo, IdxType& hi) const
  { 
    p_matrix_impl->localRowRange(lo, hi); 
  }

  /// Get the total number of rows in this matrix (specialized)
  IdxType p_rows(void) const
  { 
    return p_matrix_impl->rows(); 
  }

  /// Get the number of local rows in this matirx (specialized)
  IdxType p_localRows(void) const
  { 
    return p_matrix_impl->localRows(); 
  }

  /// Get the number of columns in this matrix (specialized)
  IdxType p_cols(void) const
  { 
    return p_matrix_impl->cols(); 
  }

  /// Get the number of local rows in this matirx (specialized)
  IdxType p_localCols(void) const
  { 
    return p_matrix_impl->localCols(); 
  }

  /// Set an individual element
  void p_setElement(const IdxType& i, const IdxType& j, const TheType& x)
  { 
    p_matrix_impl->setElement(i, j, x); 
  }

  /// Set an several element
  void p_setElements(const IdxType& n, 
                     const IdxType *i, const IdxType *j, 
                     const TheType *x)
  { 
    p_matrix_impl->setElements(n, i, j, x); 
  }

  /// Add to  an individual element
  void p_addElement(const IdxType& i, const IdxType& j, const TheType& x)
  { 
    p_matrix_impl->addElement(i, j, x); 
  }

  /// Add to  an several element
  void p_addElements(const IdxType& n, 
                     const IdxType *i, const IdxType *j, 
                     const TheType *x)
  { 
    p_matrix_impl->addElements(n, i, j, x); 
  }

  /// Get an individual element
  void p_getElement(const IdxType& i, const IdxType& j, TheType& x) const
  { 
    p_matrix_impl->getElement(i, j, x); 
  }

  /// Get an several element
  void p_getElements(const IdxType& n, 
                     const IdxType *i, const IdxType *j, 
                     TheType *x) const
  { 
    p_matrix_impl->getElements(n, i, j, x); 
  }

  /// Scale this entire MatrixT by the given value (specialized)
  void p_scale(const TheType& x)
  {
    p_matrix_impl->scale(x);
  }

  /// Shift the diagonal of this matrix by the specified value (specialized)
  void p_addDiagonal(const TheType& x)
  {
    p_matrix_impl->addDiagonal(x);
  }

  /// Make this matrix the identity matrix (specialized)
  void p_identity(void) 
  {
    p_matrix_impl->identity();
  }

  ///  Get a row and put it in a local array (specialized)
  void p_getRow(const IdxType& row, TheType *x) const
  {
    p_matrix_impl->getRow(row, x);
  }

  /// Get some rows and put them in a local array (specialized)
  void p_getRowBlock(const IdxType& nrow, const IdxType *rows, TheType *x) const
  {
    p_matrix_impl->getRowBlock(nrow, rows, x);
  }

  /// Replace all elements with their real parts
  void p_real(void)
  { 
    p_matrix_impl->real(); 
  }

  /// Replace all elements with their imaginary parts
  void p_imaginary(void)
  { 
    p_matrix_impl->imaginary(); 
  }

  /// Replace all elements with their complex gradient
  void p_conjugate(void)
  { 
    p_matrix_impl->conjugate(); 
  }

  /// Compute the matrix L<sup>2</sup> norm (specialized)
  double p_norm2(void) const
  { 
    return p_matrix_impl->norm2(); 
  }

  /// Zero all entries in the matrix (specialized)
  void p_zero(void)
  {
    p_matrix_impl->zero();
  }

  /// Make this instance ready to use
  void p_ready(void)
  { 
    p_matrix_impl->ready(); 
  }

  /// Print to named file or standard output
  void p_print(const char* filename = NULL) const
  {
    p_matrix_impl->print(filename);
  }

   /// Save, in MatLAB format, to named file (collective)
  void p_save(const char *filename) const
  {
    p_matrix_impl->save(filename);
  }

  /// Load from a named file of whatever binary format the math library uses
  void p_loadBinary(const char *filename)
  {
    p_matrix_impl->loadBinary(filename);
  }

  /// Save to named file in whatever binary format the math library uses
  void p_saveBinary(const char *filename) const
  {
    p_matrix_impl->saveBinary(filename);
  }

  /// Allow visits by implementation visitors
  void p_accept(ImplementationVisitor& visitor)
  { 
    p_matrix_impl->accept(visitor); 
  }

  /// Allow visits by implementation visitors
  void p_accept(ConstImplementationVisitor& visitor) const
  { 
    p_matrix_impl->accept(visitor); 
  }

  /// Is a MatrixT compatible with this one (throws if not)
  void p_check_compatible(const MatrixT& A) const
  {
    // FIXME: should be able to check this:
    // if (this->communicator() != A.communicator()) {
    //   throw gridpack::Exception("incompatible: communicators do not match");
    // }

    if ((this->rows() != A.rows()) || (this->cols() != A.cols())) {
      throw gridpack::Exception("incompatible: sizes do not match");
    }
    return;
  }

};

typedef MatrixT<ComplexType> ComplexMatrix;
typedef MatrixT<RealType> RealMatrix;

/// The main matrix type
typedef ComplexMatrix Matrix;


// -------------------------------------------------------------
// Matrix Operations
//
// put results in existing instance and throw when arguments are
// inconsistent
// -------------------------------------------------------------

/// Add two Matrix instances and put the result in a third
template <typename T, typename I>
void add(const MatrixT<T, I>& A, const MatrixT<T, I>& B, MatrixT<T, I>& result)
{
  result.equate(A);
  result.add(B);
}

/// Make the transpose of a Matrix and put it in another
/** 
 * 
 * 
 * @param A 
 * @param result 
 */
template <typename T, typename I>
void 
transpose(const MatrixT<T, I>& A, MatrixT<T, I>& result);

/// Get a column from the Matrix and put in specified Vector
/** 
 * 
 * 
 * @param A 
 * @param cidx 
 * @param x 
 */
template <typename T, typename I>
void 
column(const MatrixT<T, I>& A, const int& cidx, VectorT<T, I>& x);

/// Get the diagonal from a Matrix and put it in specified Vector
/** 
 * 
 * @param A 
 * @param x 
 */
template <typename T, typename I>
void 
diagonal(const MatrixT<T, I>& A, VectorT<T, I>& x);

/// Multiply two Matrix instances and put result in existing Matrix
/** 
 * 
 * 
 * @param A 
 * @param B 
 * @param result 
 */
template <typename T, typename I>
void 
multiply(const MatrixT<T, I>& A, const MatrixT<T, I>& B, MatrixT<T, I>& result);

/// Multiply a Matrix by a Vector and put result in existing Vector
/** 
 * @c A, @c x, and @c result must all have the same \ref
 * parallel::Communicator "communicator". @c x and @c result must be
 * the same size. The length of @c x must be the number of columns in
 * @c A. If these conditions are not met, an exception is thrown.
 * 
 * @param A 
 * @param x 
 * @param result 
 */
template <typename T, typename I>
void 
multiply(const MatrixT<T, I>& A, const VectorT<T, I>& x, VectorT<T, I>& result);

/// Multiply the transpose of a Matrix by a Vector and put the result in existing Vector
/** 
 * 
 * 
 * @param A 
 * @param x 
 * @param result same size and distribution as @c x
 */
template <typename T, typename I>
void 
transposeMultiply(const MatrixT<T, I>& A, const VectorT<T, I>& x, VectorT<T, I>& result);

// -------------------------------------------------------------
// Matrix Operations 
//
// these allocate new instances and throw when arguments are
// inconsistent
// -------------------------------------------------------------

/// Add two Matrix instances
/** 
 * @e Collective.
 *
 * @c A and @c B must have the same communicator and the same
 * size. Different parallel distributions and nonzero patterns are OK,
 * though.
 * 
 * @param A 
 * @param B 
 * 
 * @return pointer to new Matrix containing A+B
 */
template <typename T, typename I>
MatrixT<T, I> *add(const MatrixT<T, I>& A, const MatrixT<T, I>& B)
{
  MatrixT<T, I> *result = A.clone();
  add(A, B, *result);
  return result;
}

/// Make the transpose of a Matrix
/** 
 * @e Collective.
 *
 * 
 * 
 * @param A 
 * 
 * @return pointer to a new Matrix containing A<sup>T</sup>
 */
template <typename T, typename I>
MatrixT<T, I> *transpose(const MatrixT<T, I>& A);

/// Get a column from the Matrix and put in new Vector
/** 
 * 
 * 
 * @param A 
 * @param cidx 
 * 
 * @return 
 */
template <typename T, typename I>
VectorT<T, I> *column(const MatrixT<T, I>& A, const int& cidx)
{
  VectorT<T, I> *colv(new VectorT<T, I>(A.communicator(), A.localRows()));
  column<T, I>(A, cidx, *colv);
  return colv;
}


/// Get the diagonal from a Matrix and put in new Vector
/** 
 * @e Collective.
 * 
 * @param A 
 * 
 * @return pointer to new, allocated Vector containing the diagonal elements of @c A
 */
template <typename T, typename I>
VectorT<T, I> *diagonal(const MatrixT<T, I>& A)
{
  VectorT<T, I> *colv(new VectorT<T, I>(A.communicator(), A.localRows()));
  diagonal<T, I>(A, *colv);
  return colv;
}

/// Make a diagonal Matrix from a Vector
/** 
 * @e Collective.
 *
 * This 
 * 
 * @param x 
 * 
 * @return 
 */
template <typename T, typename I>
MatrixT<T, I> *diagonal(const VectorT<T, I>& x, 
                        const MatrixStorageType& stype = Sparse);

/// Multiply two Matrix instances and make a new one
/** 
 * 
 * 
 * @param A 
 * @param B 
 * 
 * @return 
 */
template <typename T, typename I>
MatrixT<T, I> *multiply(const MatrixT<T, I>& A, const MatrixT<T, I>& B);

/// Multiply a Matrix by a Vector and make a new Vector for the result
/** 
 * 
 * @param A 
 * @param x 
 * 
 * @return 
 */
template <typename T, typename I>
VectorT<T, I> *multiply(const MatrixT<T, I>& A, const VectorT<T, I>& x)
{
  VectorT<T, I> *result(new VectorT<T, I>(x.communicator(), A.localRows()));
  multiply<T, I>(A, x, *result);
  return result;
}

/// Multiply the transpose of a Matrix by a Vector and make a new Vector for the result
/** 
 * 
 * 
 * @param A 
 * @param x 
 * 
 * @return 
 */
template <typename T, typename I>
VectorT<T, I> *transposeMultiply(const MatrixT<T, I>& A, const VectorT<T, I>& x)
{
  VectorT<T, I> *result(new VectorT<T, I>(x.communicator(), x.localSize()));
  transposeMultiply<T, I>(A, x, *result);
  return result;
}

/// Make an identity matrix with the same ownership as the specified matrix
template <typename T, typename I>
MatrixT<T, I> *identity(const MatrixT<T, I>& A)
{
  MatrixT<T, I> *result(A.clone());
  result->identity();
  return result;
}


/// Create a new matrix containing the real part of the specified matrix
template <typename T, typename I>
MatrixT<T, I> *real(const MatrixT<T, I>& A)
{
  MatrixT<T, I> *result = A.clone();
  result->real();
  return result;
}

/// Create a new matrix containing the imaginary part of the specified matrix
template <typename T, typename I>
MatrixT<T, I> *imaginary(const MatrixT<T, I>& A)
{
  MatrixT<T, I> *result = A.clone();
  result->imaginary();
  return result;
}

/// Create a new matrix containing the complex conjugate of the specified matrix
template <typename T, typename I>
MatrixT<T, I> *conjugate(const MatrixT<T, I>& A)
{
  MatrixT<T, I> *result = A.clone();
  result->conjugate();
  return result;
}


/// Create a copy of a Matrix, possibly with a different storage type
template <typename T, typename I>
extern MatrixT<T, I> *
storageType(const MatrixT<T, I>& A, const MatrixStorageType& new_type);

/// Create a new matrix and load its contents from the specified (binary) file
template <typename T, typename I>
extern MatrixT<T, I> *
matrixLoadBinary(const parallel::Communicator&comm, const char *filename);


} // namespace math
} // namespace gridpack



#endif
