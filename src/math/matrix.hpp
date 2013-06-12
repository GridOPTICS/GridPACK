// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   matrix.hpp
 * @author William A. Perkins
 * @date   2013-06-12 10:27:11 d3g096
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
  : public parallel::Distributed,
    private utility::Uncopyable
{
public:

  /// The types of matrices that can be created
  /**
   * The gridpack::math library provides two storage schemes for
   * matrices. This is used by Matrix and MatrixImplementation
   * subclasses.
   * 
   */
  enum StorageType { Dense, Sparse };

  /// Constructor.
  Matrix(const parallel::Communicator& dist,
         const int& local_rows,
         const int& cols,
         const StorageType& storage_type = Sparse);

  /// Destructor
  ~Matrix(void);

  /// Get the total number of rows in this matrix
  int rows(void) const
  {
    return p_matrix_impl->rows();
  }

  /// Get the number of local rows in this matirx
  int local_rows(void) const
  {
    return p_matrix_impl->local_rows();
  }

  /// Get the range of global row indexes owned by this process
  void local_row_range(int& lo, int& hi) const
  {
    p_matrix_impl->local_row_range(lo, hi);
  }

  /// Get the number of columns in this matrix
  int cols(void) const
  {
    return p_matrix_impl->cols();
  }

  // /// Set an individual element
  void set_element(const int& i, const int& j, const ComplexType& x)
  {
    p_matrix_impl->set_element(i, j, x);
  }

  /// Set an several elements
  void set_elements(const int& n, const int *i, const int *j, const ComplexType *x)
  {
    p_matrix_impl->set_elements(n, i, j, x);
  }

  // /// Set all elements in a row
  // void set_row(const int& nj, const int& i, const int *j, const ComplexType *x)
  // {
  //   p_matrix_impl->set_row(nj, i, j, x);
  // }

  // /// Set all elements in a row
  // void set_region(const int& ni, const int& nj, 
  //                 const int *i, const int *j, const ComplexType *x)
  // {
  //   p_matrix_impl->set_row(ni, nj, i, j, x);
  // }

  /// Add to an individual element
  void add_element(const int& i, const int& j, const ComplexType& x)
  {
    p_matrix_impl->add_element(i, j, x);
  }

  /// Add to an several elements
  void add_elements(const int& n, const int *i, const int *j, const ComplexType *x)
  {
    p_matrix_impl->add_elements(n, i, j, x);
  }

  // /// Add to all elements in a row
  // void add_row(const int& nj, const int& i, const int *j, const ComplexType *x)
  // {
  //   p_matrix_impl->add_row(nj, i, j, x);
  // }

  /// Get an individual element
  void get_element(const int& i, const int& j, ComplexType& x) const
  {
    p_matrix_impl->get_element(i, j, x);
  }

  /// Get an several elements
  void get_elements(const int& n, const int *i, const int *j, ComplexType *x) const
  {
    p_matrix_impl->get_elements(n, i, j, x);
  }

  // /// Get all elements in a row
  // void get_row(const int& nj, const int& i, const int *j, ComplexType *x)  const
  // {
  //   p_matrix_impl->get_row(nj, i, j, x);
  // }

  // /// Get all elements in a row
  // void get_region(const int& ni, const int& nj, 
  //                 const int *i, const int *j, ComplexType *x) const
  // {
  //   p_matrix_impl->get_row(ni, nj, i, j, x);
  // }

  /// Indicate the matrix is ready to use
  void ready(void)
  {
    p_matrix_impl->ready();
  }

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

  /// Make an exact replica of this instance
  Matrix *clone(void) const
  {
    MatrixImplementation *pimpl_clone =
      this->p_matrix_impl->clone();
    Matrix *result = new Matrix(pimpl_clone);
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
  void equate(const Matrix& A);

  /// Scale the entire Matrix by the given value
  void scale(const ComplexType& x);

  /// Multiply the diagonal by the specified vector
  void multiply_diagonal(const Vector& x);

  /// Add another matrix to this one, in place
  void add(const Matrix& A);

  /// Make this matrix the identity matrix
  void identity(void);

  /// Zero all entries in the matrix
  void zero(void);

  // -------------------------------------------------------------
  // Matrix Operations 
  //
  // all allocate new instances and throw when arguments are
  // inconsistent
  // -------------------------------------------------------------
  friend Matrix *factorize(const Matrix& A);
  friend Matrix *inverse(const Matrix& A);
  // friend Matrix *reorder(const Matrix& A, const Reordering& r);
  friend Matrix *identity(const Matrix& A);

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
extern Matrix *add(const Matrix& A, const Matrix& B);

/// Make the transpose of a Matrix
extern Matrix *transpose(const Matrix& A);

/// Get a column from the Matrix and put in new Vector
extern Vector *column(const Matrix& A, const int& cidx);

/// Get the diagonal from from a Matrix and put in new Vector
extern Vector *diagonal(const Matrix& A);

/// Multiply two Matrix instances and make a new one
extern Matrix *multiply(const Matrix& A, const Matrix& B);

/// Multiply a Matrix by a Vector and make a new Vector for the result
extern Vector *multiply(const Matrix& A, const Vector& x);

// -------------------------------------------------------------
// Matrix Operations
//
// put results in existing instance and throw when arguments are
// inconsistent
// -------------------------------------------------------------

/// Add two Matrix instances and put the result in a third
extern void add(const Matrix& A, const Matrix& B, Matrix& result);

/// Make the transpose of a Matrix and put it in another
extern void transpose(const Matrix& A, Matrix& result);

/// Get a column from the Matrix and put in specified Vector
extern void column(const Matrix& A, const int& cidx, Vector& x);

/// Get the diagonal from a Matrix and put it in specified Vector
extern void diagonal(const Matrix& A, Vector& x);

/// Multiply two Matrix instances and put result in existing Matrix
extern void multiply(const Matrix& A, const Matrix& B, Matrix& result);

/// Multiply a Matrix by a Vector and put result in existing Vector
extern void multiply(const Matrix& A, const Vector& x, Vector& result);


} // namespace math
} // namespace gridpack



#endif
