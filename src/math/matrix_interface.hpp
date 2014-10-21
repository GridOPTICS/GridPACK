// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   matrix_interface.hpp
 * @author William A. Perkins
 * @date   2014-10-21 13:30:33 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 21, 2014 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------

#ifndef _matrix_interface_hpp_
#define _matrix_interface_hpp_

namespace gridpack {
namespace math {


// -------------------------------------------------------------
//  class BaseMatrixInterface
// -------------------------------------------------------------
template<typename T, typename I = int>
class BaseMatrixInterface 
{
public:

  typedef T TheType;            /**< The numeric type used. */
  typedef I IdxType;            /**< The size/index type used. */

  /// Default constructor.
  BaseMatrixInterface(void)
  {}

  /// Destructor
  ~BaseMatrixInterface(void)
  {}

  /// Get the global index range of the locally owned rows
  /** 
   * @e Local.
   * 
   * 
   * 
   * @param lo first (0-based) index of locally owned rows
   * @param hi one more than the last (0-based) index of locally owned rows
   */
  void localRowRange(IdxType& lo, IdxType& hi) const 
  {
    this->p_localRowRange(lo, hi);
  }

  /// Get the total number of rows in this matrix
  /** 
   * @e Local.
   * 
   * 
   * @return total number of rows
   */
  IdxType rows(void) const
  {
    return this->p_rows();
  }

  /// Get the number of local rows in this matirx
  IdxType localRows(void) const
  /** 
   * @e Local.
   * 
   * 
   * @return number of rows owned by this process
   */
  {
    return this->p_localRows();
  }

  /// Get the number of columns in this matrix
  /** 
   * @e Local.
   * 
   * 
   * @return number of columns in matrix
   */
  IdxType cols(void) const
  {
    return this->p_cols();
  }

  /// Get the number of local columns in this matirx
  IdxType localCols(void) const
  {
    return this->p_localCols();
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
  void setElement(const IdxType& i, const IdxType& j, const TheType& x)
  {
    this->p_setElement(i, j, x);
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
  void setElements(const IdxType& n, const IdxType *i, const IdxType *j, const TheType *x)
  {
    this->p_setElements(n, i, j, x);
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
  void addElement(const IdxType& i, const IdxType& j, const TheType& x)
  {
    this->p_addElement(i, j, x);
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
  void addElements(const IdxType& n, const IdxType *i, const IdxType *j, const TheType *x)
  {
    this->p_addElements(n, i, j, x);
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
  void getElement(const IdxType& i, const IdxType& j, TheType& x) const
  {
    this->p_getElement(i, j, x);
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
  void getElements(const IdxType& n, const IdxType *i, const IdxType *j, TheType *x) const
  {
    this->p_getElements(n, i, j, x);
  }

  /// Replace all elements with their real parts
  void real(void)
  {
    this->p_real();
  }

  /// Replace all elements with their imaginary parts
  void imaginary(void)
  {
    this->p_imaginary();
  }

  /// Replace all elements with their complex gradient
  void conjugate(void)
  {
    this->p_conjugate();
  }

  /// Compute the matrix L<sup>2</sup> norm
  /** 
   * @e Collective.
   *
   * The vector L<sup>2</sup>, or Euclidian, norm is computed as
   * \f[
   *   \left\| \mathbf{A} \right\| ~ = ~ \sqrt{\sum_{ij} A_{ij}^{2}}
   * \f]
   * 
   * @return  L<sup>2</sup> norm of this matrix
   */
  double norm2(void) const
  {
    return this->p_norm2();
  }

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
    this->p_ready();
  }

protected:

  /// Get the global index range of the locally owned rows (specialized)
  virtual void p_localRowRange(IdxType& lo, IdxType& hi) const = 0;

  /// Get the total number of rows in this matrix (specialized)
  virtual IdxType p_rows(void) const = 0;

  /// Get the number of local rows in this matirx (specialized)
  virtual IdxType p_localRows(void) const = 0;

  /// Get the number of columns in this matrix (specialized)
  virtual IdxType p_cols(void) const = 0;

  /// Get the number of local rows in this matirx (specialized)
  virtual IdxType p_localCols(void) const = 0;

  /// Set an individual element
  virtual void p_setElement(const IdxType& i, const IdxType& j, const TheType& x) = 0;

  /// Set an several element
  virtual void p_setElements(const IdxType& n, const IdxType *i, const IdxType *j, 
                             const TheType *x) = 0;

  /// Add to  an individual element
  virtual void p_addElement(const IdxType& i, const IdxType& j, const TheType& x) = 0;

  /// Add to  an several element
  virtual void p_addElements(const IdxType& n, const IdxType *i, const IdxType *j, 
                             const TheType *x) = 0;

  /// Get an individual element
  virtual void p_getElement(const IdxType& i, const IdxType& j, TheType& x) const = 0;

  /// Get an several element
  virtual void p_getElements(const IdxType& n, const IdxType *i, const IdxType *j, 
                             TheType *x) const = 0;

  /// Replace all elements with their real parts
  virtual void p_real(void) = 0;

  /// Replace all elements with their imaginary parts
  virtual void p_imaginary(void) = 0;

  /// Replace all elements with their complex gradient
  virtual void p_conjugate(void) = 0;

  /// Compute the matrix L<sup>2</sup> norm (specialized)
  virtual double p_norm2(void) const = 0;

  /// Make this instance ready to use
  virtual void p_ready(void) = 0;

};


} // namespace math
} // namespace gridpack

#endif
