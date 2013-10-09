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
 * @date   2013-10-09 12:23:55 d3g096
 * 
 * @brief  
 * 
 * 
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
class Matrix 
  : public parallel::WrappedDistributed,
    private utility::Uncopyable
{
public:

  /// The types of matrices that can be created
  /**
   * The gridpack::math library provides two storage schemes for
   * matrices. This is used by Matrix and MatrixImplementation
   * subclasses.
   *
   * The actual storage scheme and memory used is dependent upon the
   * underlying math library implementation.
   * 
   */
  enum StorageType { 
    Dense,                      /**< dense matrix storage scheme */
    Sparse                      /**< sparse matrix storage scheme */
  };

  /// Constructor.
  /** 
   * A Matrix must be instantiated simulutaneously on all processes
   * involved in the specified \ref parallel::Communicator
   * "communicator". Each process in the communicator will own the
   * number of rows requested.
   * 
   * @param dist parallel environment
   * @param local_rows matrix rows to be owned by the local process
   * @param cols total number of columns (same on all processes)
   * @param storage_type specify dense or sparse storage 
   * 
   * @return 
   */
  Matrix(const parallel::Communicator& dist,
         const int& local_rows,
         const int& cols,
         const StorageType& storage_type = Sparse);

  /// Destructor
  /** 
   * A matrix must be destroyed simulutaneously on all processes
   * involved in the \ref parallel::Communicator "communicator" used
   * to instantiate it.
   */
  ~Matrix(void);

  /// Get the global number of rows in this matrix
  /** 
   * @e Local.
   * 
   * 
   * @return total number of rows
   */
  int rows(void) const
  {
    return p_matrix_impl->rows();
  }

  /// Get the number of local rows in this matirx
  /** 
   * @e Local.
   * 
   * 
   * @return number of rows owned by this process
   */
  int localRows(void) const
  {
    return p_matrix_impl->localRows();
  }

  /** 
   * @deprecated does not meet coding standards, use localRows()
   * 
   * 
   * @return 
   */
  int local_rows(void) const
  {
    return this->localRows();
  }

  /// Get the range of global row indexes owned by this process
  /** 
   * @e Local.
   * 
   * 
   * 
   * @param lo first (0-based) index of locally owned rows
   * @param hi one more than the last (0-based) index of locally owned rows
   */
  void localRowRange(int& lo, int& hi) const
  {
    p_matrix_impl->localRowRange(lo, hi);
  }

  /** 
   * @deprecated does not meet coding standards, use localRowRange()
   * 
   * 
   * @param lo 
   * @param hi 
   */
  void local_row_range(int& lo, int& hi) const
  {
    this->localRowRange(lo, hi);
  }

  /// Get the global number of columns in this matrix
  /** 
   * @e Local.
   * 
   * 
   * @return number of columns in matrix
   */
  int cols(void) const
  {
    return p_matrix_impl->cols();
  }

  /// Set an individual element
  /** 
   * @e Local.
   *
   * This overwrites the value at the specified @c i row and @c j
   * column.  ready() must be called after all setElement() calls and
   * before using the matrix.
   * 
   * @param i global (0-based) row index
   * @param j global (0-based) column index
   * @param x value to place in matrix
   */
  void setElement(const int& i, const int& j, const ComplexType& x)
  {
    p_matrix_impl->setElement(i, j, x);
  }

  /** 
   * @deprecated does not meet coding standards, use setElement()
   * 
   * 
   * @param i 
   * @param j 
   * @param x 
   */
  void set_element(const int& i, const int& j, const ComplexType& x)
  {
    p_matrix_impl->setElement(i, j, x);
  }

  /// Set an several elements
  /** 
   * @e Local.
   *
   * This overwrites values at several locations in the
   * matrix. ready() must be called after all setElements() calls and
   * before using the matrix.
   * 
   * @param n number of values to place in 
   * @param i array of @c n global, 0-based row indexes
   * @param j array of @c n global, 0-based column indexes
   * @param x array of @c n values to replace existing matrix elements
   */
  void setElements(const int& n, const int *i, const int *j, const ComplexType *x)
  {
    p_matrix_impl->setElements(n, i, j, x);
  }

  /// Add to an individual element
  /** 
   * @e Local.
   * 
   * 
   * @param i global, 0-based row index
   * @param j global, 0-based column index
   * @param x value to add to matrix element
   */
  void addElement(const int& i, const int& j, const ComplexType& x)
  {
    p_matrix_impl->addElement(i, j, x);
  }

  /** 
   * @deprecated does not meet coding standards, use addElement()
   * 
   * @param i 
   * @param j 
   * @param x 
   */
  void add_element(const int& i, const int& j, const ComplexType& x)
  {
    this->addElement(i, j, x);
  }

  /// Add to an several elements
  /** 
   * @e Local.
   * 
   * @param n number of elements to update
   * @param i array of @c n global, 0-based row indexes
   * @param j array of @c n global, 0-based column indexes
   * @param x array of @c n values to add to existing matrix elements
   */
  void addElements(const int& n, const int *i, const int *j, const ComplexType *x)
  {
    p_matrix_impl->addElements(n, i, j, x);
  }

  /** 
   * @deprecated does not meet coding standards, use addElements()
   * 
   * @param n 
   * @param i 
   * @param j 
   * @param x 
   */
  void add_elements(const int& n, const int *i, const int *j, const ComplexType *x)
  {
    this->addElements(n, i, j, x);
  }

  /// Get an individual element
  /** 
   * @c Local.
   * 
   * Only local elements may be retrieved using this method.
   * 
   * @param i global, 0-based row index
   * @param j global, 0-based column index
   * @param x variable in which to place retrieved value
   */
  void getElement(const int& i, const int& j, ComplexType& x) const
  {
    p_matrix_impl->getElement(i, j, x);
  }

  /** 
   * @deprecated does not meet coding standards, use getElement
   * 
   * 
   * @param i 
   * @param j 
   * @param x 
   */
  void get_element(const int& i, const int& j, ComplexType& x) const
  {
    this->getElement(i, j, x);
  }

  /// Get an several elements
  /** 
   * @c Local.
   *
   * Only local elements may be retrieved using this method.
   * 
   * @param n number of elements to retrieve
   * @param i array of @c n global, 0-based row indexes
   * @param j array of @c n global, 0-based column indexes
   * @param x array of @c n values in which retrieved values are to be placed
   */
  void getElements(const int& n, const int *i, const int *j, ComplexType *x) const
  {
    p_matrix_impl->getElements(n, i, j, x);
  }

  /** 
   * @deprecated does not meet coding standards, use getElements
   * 
   * @param n 
   * @param i 
   * @param j 
   * @param x 
   */
  void get_elements(const int& n, const int *i, const int *j, ComplexType *x) const
  {
    this->getElements(n, i, j, x);
  }

  /// Make this matrix the identity matrix
  /** 
   * @e Collective
   *
   * 
   * 
   */
  void identity(void);

  /// Indicate the matrix is ready to use
  /** 
   * @e Collective.
   *
   * This is used to indicate that the matrix is ready to use.  This
   * must be called after @e all setElement() or addElement() calls
   * and before the vector is used for any operation.
   */
  void ready(void)
  {
    p_matrix_impl->ready();
  }

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
  Matrix *clone(void) const
  {
    MatrixImplementation *pimpl_clone =
      this->p_matrix_impl->clone();
    Matrix *result = new Matrix(pimpl_clone);
    return result;
  }

  /// Print to named file or standard output
  /** 
   * @e Collective
   * 
   * @param filename name of file to write, or NULL for standard output
   */
  void print(const char* filename = NULL) const;

  /// Save, in MatLAB format, to named file
  /** 
   * @e Collective
   * 
   * @param filename name of file to write
   */
  void save(const char *filename) const;

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
  void equate(const Matrix& A);

  /// Scale the entire Matrix by the given value
  /** 
   * @e Collective.
   * 
   * @param x factor by which all elements in the matrix are multiplied
   */
  void scale(const ComplexType& x);

  /// Multiply the diagonal by the specified vector
  /** 
   * @e Collective.
   * 
   * @param x factor by which all diagonal elements in the matrix are multiplied
   */
  void multiplyDiagonal(const Vector& x);

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
  void add(const Matrix& A);

  /// Zero all entries in the matrix
  /** 
   * @e Collective.
   * 
   */
  void zero(void);

  // -------------------------------------------------------------
  // Matrix Operations 
  //
  // all allocate new instances and throw when arguments are
  // inconsistent
  // -------------------------------------------------------------
  // friend Matrix *factorize(const Matrix& A);
  // friend Matrix *inverse(const Matrix& A);
  // friend Matrix *reorder(const Matrix& A, const Reordering& r);
  // friend Matrix *identity(const Matrix& A);
  friend Matrix *multiply(const Matrix& A, const Matrix& B);

protected:

  /// The actual implementation
  boost::scoped_ptr<MatrixImplementation> p_matrix_impl;

  /// Construct with an existing (allocated) implementation 
  Matrix(MatrixImplementation *impl);

  /// Is a Matrix compatible with this one (throws if not)
  void p_check_compatible(const Matrix& A) const;

};

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
extern Matrix *add(const Matrix& A, const Matrix& B);

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
extern Matrix *transpose(const Matrix& A);

/// Get a column from the Matrix and put in new Vector
/** 
 * 
 * 
 * @param A 
 * @param cidx 
 * 
 * @return 
 */
extern Vector *column(const Matrix& A, const int& cidx);

/// Get the diagonal from from a Matrix and put in new Vector
/** 
 * @e Collective.
 * 
 * @param A 
 * 
 * @return 
 */
extern Vector *diagonal(const Matrix& A);

/// Multiply two Matrix instances and make a new one
/** 
 * 
 * 
 * @param A 
 * @param B 
 * 
 * @return 
 */
extern Matrix *multiply(const Matrix& A, const Matrix& B);

/// Multiply a Matrix by a Vector and make a new Vector for the result
/** 
 * 
 * @param A 
 * @param x 
 * 
 * @return 
 */
extern Vector *multiply(const Matrix& A, const Vector& x);

/// Make an identity matrix with the same ownership as the specified matrix
extern Matrix *identity(const Matrix& A);

// -------------------------------------------------------------
// Matrix Operations
//
// put results in existing instance and throw when arguments are
// inconsistent
// -------------------------------------------------------------

/// Add two Matrix instances and put the result in a third
extern void add(const Matrix& A, const Matrix& B, Matrix& result);

/// Make the transpose of a Matrix and put it in another
/** 
 * 
 * 
 * @param A 
 * @param result 
 */
extern void transpose(const Matrix& A, Matrix& result);

/// Get a column from the Matrix and put in specified Vector
/** 
 * 
 * 
 * @param A 
 * @param cidx 
 * @param x 
 */
extern void column(const Matrix& A, const int& cidx, Vector& x);

/// Get the diagonal from a Matrix and put it in specified Vector
/** 
 * 
 * @param A 
 * @param x 
 */
extern void diagonal(const Matrix& A, Vector& x);

/// Multiply two Matrix instances and put result in existing Matrix
/** 
 * 
 * 
 * @param A 
 * @param B 
 * @param result 
 */
extern void multiply(const Matrix& A, const Matrix& B, Matrix& result);

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
extern void multiply(const Matrix& A, const Vector& x, Vector& result);

} // namespace math
} // namespace gridpack



#endif
