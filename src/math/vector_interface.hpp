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
 * @file   vector_interface.hpp
 * @author William A. Perkins
 * @date   2014-10-21 11:19:20 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 20, 2014 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _vector_interface_hpp_
#define _vector_interface_hpp_

namespace gridpack {
namespace math {

// -------------------------------------------------------------
// Abstract class VectorInterface
// -------------------------------------------------------------
/**
 * This defines the basic interface for vectors.  
 * 
 */
template<typename T, typename I = int> 
class BaseVectorInterface 
{
public:

  typedef T TheType;            /**< The numeric type used */
  typedef I IdxType;            /**< The size/index type used */

  /// Default constructor.
  BaseVectorInterface(void)
  {}

  /// Destructor
  ~BaseVectorInterface(void)
  {}

  /// Get the global length
  /** 
   * 
   * @return global vector length
   */
  IdxType size(void) const
  {
    return this->p_size();
  }

  /// Get the local length
  /** 
   * 
   * @return local vector length
   */
  IdxType localSize(void) const
  {
    return this->p_localSize();
  }

  /// Get the local min/max global indexes
  /** 
   * @e Local.
   * 
   * The minimum index in the first global (0-based) index owned by
   * the local process. The maximum is one more than the last global
   * index owned. An example of usage:
   * 
     \code{.cpp}
     int lo, hi;
     my_vector.local_index_range(lo, hi);
     for (int i = lo; i < hi; ++i) {
       ComplexType x;
       x = ...;
       v.setElement(i, x);
     }
     \endcode
   * 
   * @param lo first (0-based) index of locally owned elements
   * @param hi one more than the last (0-based) index of locally owned elements
   */
  void localIndexRange(IdxType& lo, IdxType& hi) const
  {
    return this->p_localIndexRange(lo, hi);
  }

  /// Set an individual element
  /** 
   * @e Local.
   * 
   * This overwrites the value at the specified index.  ready() must
   * be called after all setElement() calls and before using the
   * vector.
   * 
   * @param i element global (0-based) index 
   * @param x value to place in vector
   */
  void setElement(const IdxType& i, const TheType& x)
  {
    this->p_setElement(i, x);
  }

  /// Set a range of elements (lo to hi-1)
  /** 
   * @e Local.
   * 
   * An example that fills the locally owned part of the vector:
   * \code{.cpp}
   * Vector v(...);
   * int lo, hi;
   * v.local_index_range(lo, hi);
   * std::vector<ComplexType> x(v.local_size())
   * // fill x with appropriate values
   * v.setElement_range(lo, hi, &x[0]);
   * \endcode
   * 
   * @param lo lowest global (0-based) index to fill
   * @param hi one more than the highest global (0-based) index to fill
   * @param x array of hi - lo values
   */
  void setElementRange(const IdxType& lo, const IdxType& hi, TheType *x)
  {
    this->p_setElementRange(lo, hi, x);
  }

  /// Set an several elements
  /** 
   * @e Local.
   * 
   * This places (overwrites) several elements, with arbitrary
   * indexes, in the vector.  ready() must be called after all
   * setElement() calls and before using the vector.
   * 
   * @param n number of elements to place in vector
   * @param i pointer to an array of @c n global (0-based) indexes
   * @param x pointer to an arry 
   */
  void setElements(const IdxType& n, const IdxType *i, const TheType *x)
  {
    this->p_setElements(n, i, x);
  }

  /// Add to an individual element
  /** 
   * @e Local.
   * 
   * 
   * 
   * @param i 
   * @param x 
   */
  void addElement(const IdxType& i, const TheType& x)
  {
    this->p_addElement(i, x);
  }

  /// Add to an several elements
  /** 
   * @e Local.
   * 
   * 
   * 
   * @param n 
   * @param i 
   * @param x 
   */
  void addElements(const IdxType& n, const IdxType *i, const TheType *x)
  {
    this->p_addElements(n, i, x);
  }

  /// Get an individual element
  /** 
   * @e Local.
   * 
   * Only local values may be retreived with this method.  
   * 
   * @param i element's global, 0-based index
   * @param x element's value
   */
  void getElement(const IdxType& i, TheType& x) const
  {
    this->p_getElement(i, x);
  }

  /// Get an several elements
  /** 
   * @e Local.
   * 
   * Only local values may be retreived with this method.  
   * 
   * 
   * @param n number of elements to get
   * @param i array of @c n global, 0-based indexes
   * @param x array of @c n values
   */
  void getElements(const IdxType& n, const IdxType *i, TheType *x) const
  {
    this->p_getElements(n, i, x);
  }

  /// Get a range of elements (lo to hi-1)
  /** 
   * @e Local.
   *
   * The elements whose indexes range from @c lo to @c hi-1 are
   * retreived.  Only local values may be retreived with this method.
   * The length of the @c x array must be at least @c "hi - lo"
   * 
   * @param lo lowest element index
   * @param hi one more than the highest element index
   * @param x existing array of values to be filled with element values
   */
  void getElementRange(const IdxType& lo, const IdxType& hi, TheType *x) const
  {
    this->p_getElementRange(lo, hi, x);
  }

  /// Get all of vector elements (on all processes)
  /** 
   * @e Collective.
   *
   * This is an all gather operation and consequently will be slow.
   * All of this vector's elements are gathered on all processes and
   * placed in the in the specified array.
   * 
   * 
   * @param x array of size() length to be filled with values
   */
  void getAllElements(TheType *x) const
  {
    this->p_getAllElements(x);
  }


  /// Make all the elements zero
  void zero(void)
  {
    this->p_zero();
  }

  /// Make all the elements the specified value
  void fill(const TheType& v)
  {
    this->p_fill(v);
  }

  /// Compute the vector L1 norm (sum of absolute value)
  double norm1(void) const
  {
    return this->p_norm1();
  }

  /// Compute the vector L2 norm (root of sum of squares)
  double norm2(void) const
  {
    return this->p_norm2();
  }

  /// Compute the infinity (or maximum) norm
  double normInfinity(void) const
  {
    return this->p_normInfinity();
  }

  /// Replace all elements with its absolute value (complex magnitude) 
  void abs(void)
  {
    this->p_abs();
  }

  /// Replace all elements with their real part
  void real(void)
  {
    this->p_real();
  }

  /// Replace all elements with their imaginary part
  void imaginary(void)
  {
    this->p_imaginary();
  }

  /// Replace all elements with their complex conjugate
  void conjugate(void)
  {
    this->p_conjugate();
  }

  /// Replace all elements with its exponential
  void exp(void)
  {
    this->p_exp();
  }

  /// Make this instance ready to use
  void ready(void)
  {
    this->p_ready();
  }

protected:

  /// Get the global vector length (specialized)
  virtual IdxType p_size(void) const = 0;

  /// Get the size of the vector local part (specialized)
  virtual IdxType p_localSize(void) const = 0;

  /// Get the local min/max global indexes (specialized)
  virtual void p_localIndexRange(IdxType& lo, IdxType& hi) const = 0;

  /// Set an individual element (specialized)
  virtual void p_setElement(const IdxType& i, const TheType& x) = 0;

  /// Set an several elements (specialized)
  virtual void p_setElements(const IdxType& n, const IdxType *i, const TheType *x) = 0;

  /// Get a range of elements (lo to hi-1) (specialized)
  virtual void p_setElementRange(const IdxType& lo, const IdxType& hi, TheType *x) = 0;

  /// Add to an individual element (specialized)
  virtual void p_addElement(const IdxType& i, const TheType& x) = 0;

  /// Add to an several elements (specialized)
  virtual void p_addElements(const IdxType& n, const IdxType *i, const TheType *x) = 0;

  /// Get an individual element (specialized)
  virtual void p_getElement(const IdxType& i, TheType& x) const = 0;

  /// Get an several elements (specialized)
  virtual void p_getElements(const IdxType& n, const IdxType *i, TheType *x) const = 0;

  /// Get a range of elements (lo to hi-1) (specialized)
  virtual void p_getElementRange(const IdxType& lo, const IdxType& hi, TheType *x) const = 0;

  /// Get all of vector elements (on all processes)
  virtual void p_getAllElements(TheType *x) const = 0;

  /// Make all the elements zero (specialized)
  virtual void p_zero(void) = 0;

  /// Fill all the elements with the specified value (specialized)
  virtual void p_fill(const TheType& v) = 0;

  /// Compute the vector L1 norm (sum of absolute value) (specialized)
  virtual double p_norm1(void) const = 0;

  /// Compute the vector L2 norm (root of sum of squares) (specialized)
  virtual double p_norm2(void) const = 0;

  /// Compute the vector infinity (or maximum) norm (specialized)
  virtual double p_normInfinity(void) const = 0;

  /// Replace all elements with its absolute value (specialized) 
  virtual void p_abs(void) = 0;

   /// Replace all elements with their real part (specialized)
  virtual void p_real(void) = 0;
 
   /// Replace all elements with their imaginary part (specialized)
  virtual void p_imaginary(void) = 0;
 
   /// Replace all elements with their complex conjugate
  virtual void p_conjugate(void) = 0;
 
  /// Replace all elements with its exponential (specialized)
  virtual void p_exp(void) = 0;

  /// Make this instance ready to use
  virtual void p_ready(void) = 0;
};





} // namespace math
} // namespace gridpack


#endif
