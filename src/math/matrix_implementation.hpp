// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   matrix_implementation.h
 * @author William A. Perkins
 * @date   2014-02-19 11:46:28 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 25, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _matrix_implementation_h_
#define _matrix_implementation_h_

#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/utilities/complex.hpp>

namespace gridpack {
namespace math {

class ImplementationVisitor;
class ConstImplementationVisitor;

// -------------------------------------------------------------
//  class MatrixImplementation
// -------------------------------------------------------------
class MatrixImplementation 
  : private utility::Uncopyable,
    public parallel::Distributed
{
public:

  /// Default constructor.
  MatrixImplementation(const parallel::Communicator& comm);

  /// Destructor
  virtual ~MatrixImplementation(void);

  /// Get the global index range of the locally owned rows
  void localRowRange(int& lo, int& hi) const 
  {
    this->p_localRowRange(lo, hi);
  }

  /// Get the total number of rows in this matrix
  int rows(void) const
  {
    return this->p_rows();
  }

  /// Get the number of local rows in this matirx
  int localRows(void) const
  {
    return this->p_localRows();
  }

  /// Get the number of columns in this matrix
  int cols(void) const
  {
    return this->p_cols();
  }

  /// Get the number of local columns in this matirx
  int localCols(void) const
  {
    return this->p_localCols();
  }

  /// Set an individual element
  void setElement(const int& i, const int& j, const ComplexType& x)
  {
    this->p_setElement(i, j, x);
  }

  /// Set an several elements
  void setElements(const int& n, const int *i, const int *j, const ComplexType *x)
  {
    this->p_setElements(n, i, j, x);
  }

  // /// Set all elements in a row
  // void set_row(const int& nj, const int& i, const int *j, const ComplexType *x)
  // {
  //   this->p_set_row(nj, i, j, x);
  // }

  // /// Set all elements in a row
  // void set_region(const int& ni, const int& nj, const int *i, const int *j, const ComplexType *x)
  // {
  //   this->p_set_region(ni, nj, i, j, x);
  // }

  /// Add to an individual element
  void addElement(const int& i, const int& j, const ComplexType& x)
  {
    this->p_addElement(i, j, x);
  }

  /// Add to an several elements
  void addElements(const int& n, const int *i, const int *j, const ComplexType *x)
  {
    this->p_addElements(n, i, j, x);
  }

  // /// Add to all elements in a row
  // void add_row(const int& nj, const int& i, const int *j, const ComplexType *x)
  // {
  //   this->p_add_row(nj, i, j, x);
  // }

  /// Get an individual element
  void getElement(const int& i, const int& j, ComplexType& x) const
  {
    this->p_getElement(i, j, x);
  }

  /// Get an several elements
  void getElements(const int& n, const int *i, const int *j, ComplexType *x) const
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

  // /// Get all elements in a row
  // void get_row(const int& nj, const int& i, const int *j, const ComplexType *x)
  // {
  //   this->p_get_row(nj, i, j, x);
  // }

  // /// Get all elements in a row
  // void get_region(const int& ni, const int& nj, 
  //                 const int *i, const int *j, const ComplexType *x)
  // {
  //   this->p_get_region(ni, nj, i, j, x);
  // }

  /// Compute the matrix L<sup>2</sup> norm
  ComplexType norm2(void) const
  {
    return this->p_norm2();
  }

  /// Indicate the matrix is ready to use
  void ready(void)
  {
    this->p_ready();
  }

  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor)
  {
    this->p_accept(visitor);
  }

  /// Allow visits by const implemetation visitor
  void accept(ConstImplementationVisitor& visitor) const
  {
    this->p_accept(visitor);
  }

  /// Make an exact replica of this instance
  MatrixImplementation *clone(void) const
  {
    return this->p_clone();
  }


protected:

  /// Get the global index range of the locally owned rows (specialized)
  virtual void p_localRowRange(int& lo, int& hi) const = 0;

  /// Get the total number of rows in this matrix (specialized)
  virtual int p_rows(void) const = 0;

  /// Get the number of local rows in this matirx (specialized)
  virtual int p_localRows(void) const = 0;

  /// Get the number of columns in this matrix (specialized)
  virtual int p_cols(void) const = 0;

  /// Get the number of local rows in this matirx (specialized)
  virtual int p_localCols(void) const = 0;

  /// Set an individual element
  virtual void p_setElement(const int& i, const int& j, const ComplexType& x) = 0;

  /// Set an several element
  virtual void p_setElements(const int& n, const int *i, const int *j, 
                              const ComplexType *x) = 0;

  /// Add to  an individual element
  virtual void p_addElement(const int& i, const int& j, const ComplexType& x) = 0;

  /// Add to  an several element
  virtual void p_addElements(const int& n, const int *i, const int *j, 
                              const ComplexType *x) = 0;

  // /// Add to  all elements in a row
  // virtual void p_add_row(const int& i, const int *j, const ComplexType *x) = 0;

  /// Get an individual element
  virtual void p_getElement(const int& i, const int& j, ComplexType& x) const = 0;

  /// Get an several element
  virtual void p_getElements(const int& n, const int *i, const int *j, 
                              ComplexType *x) const = 0;

  /// Replace all elements with their real parts
  virtual void p_real(void) = 0;

  /// Replace all elements with their imaginary parts
  virtual void p_imaginary(void) = 0;

  /// Replace all elements with their complex gradient
  virtual void p_conjugate(void) = 0;

  /// Compute the matrix L<sup>2</sup> norm (specialized)
  virtual ComplexType p_norm2(void) const = 0;

  /// Make this instance ready to use
  virtual void p_ready(void) = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ImplementationVisitor& visitor) = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ConstImplementationVisitor& visitor) const = 0;

  /// Make an exact replica of this instance (specialized)
  virtual MatrixImplementation *p_clone(void) const = 0;
};

} // namespace utility
} // namespace gridpack



#endif
