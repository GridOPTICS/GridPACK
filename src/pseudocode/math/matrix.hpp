// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   matrix.hpp
 * @author William A. Perkins
 * @date   Fri Mar 22 11:57:00 2013
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
#include "gridpack/parallel/Distributed.hpp"
#include "gripack/math/matrix_implementation.hpp"
#include "gripack/math/vector.hpp"


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
 * serial vector is created, otherwise it's parallel. 
 * 
 */
class Matrix 
  : public parallel::Distributed,
    private UnCopyable
{
public:

  /// The types of matrices that can be created
  /**
   * The gridpack::math library provides two storage schemes for
   * matrices. This is used by Matrix and MatrixImplementation
   * subclasses.
   * 
   */
  enum StorageType {
    DenseMatrix,
    SparseMatrix
  };

  /// Default constructor.
  Matrix(const parallel::Distribution& dist, const StorageType& storage_type = SparseMatrix);

  /// Destructor
  virtual ~Matrix(void);

  /// Get the number of rows in this matrix
  int rows(void) const
  {
    matrix_impl_->rows();
  }

  /// Get the number of columns in this matrix
  int cols(void) const
  {
    matrix_impl_->cols();
  }

  /// Set an individual element
  void set_element(const int& i, const int& j, const complex_type& x)
  {
    matrix_impl_->set_element_(i, j, x);
  }

  /// Set an several elements
  void set_elements(cont int& n, const int *i, const int *j, const complex_type *x)
  {
    matrix_impl_->set_elements_(n, i, j, x);
  }

  /// Set all elements in a row
  void set_row(const int& nj, const int& i, const int *j, const complex_type *x)
  {
    matrix_impl_->set_row_(nj, i, j, x);
  }

  /// Set all elements in a row
  void set_region(const int& ni, const int& nj, const int *i, const int *j, const complex_type *x)
  {
    matrix_impl_->set_row_(ni, nj, i, j, x);
  }

  /// Add to an individual element
  void add_element(const int& i, const int& j, const complex_type& x)
  {
    matrix_impl_->add_element_(i, j, x);
  }

  /// Add to an several elements
  void add_elements(const int& n, const int *i, const int *j, const complex_type *x)
  {
    matrix_impl_->add_elements_(n, i, j, x);
  }

  /// Add to all elements in a row
  void add_row(const int& nj, const int& i, const int *j, const complex_type *x)
  {
    matrix_impl_->add_row_(nj, i, j, x);
  }

  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor) const
  {
    matrix_impl_->accept(visitor);
  }

  // -------------------------------------------------------------
  // Matrix Operation Methods
  // -------------------------------------------------------------

  /// Get a part of the matrix (new instance allocated)
  Matrix *submatrix(const int& istart, const int& jstart,
                    const int& iend, const int& jend) const;

  // -------------------------------------------------------------
  // In-Place Matrix Operation Methods (change this instance)
  // -------------------------------------------------------------
  void scale(const complex_type& x);
  void multiply_diagonal(const Vector& x);
  void add(const Matrix& A);

  // -------------------------------------------------------------
  // Matrix Operations (all allocate new instances)
  // -------------------------------------------------------------
  friend Matrix *add(const Matrix& A, const Matrix& B);
  friend Matrix *multiply(const Matrix& A, const Matrix& B);
  friend Vector *multiply(const Matrix& A, const Vector& x);
  friend Matrix *factorize(const Matrix& A);
  friend Matrix *transpose(const Matrix& A);
  friend Matrix *inverse(const Matrix& A);
  friend Matrix *reorder(const Matrix& A, const Reordering& r);
  

protected:

  boost::scoped_ptr<MatrixImplementation> matrix_impl_;

};

} // namespace math
} // namespace gridpack



#endif
