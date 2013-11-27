// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   vector_implementation.h
 * @author William A. Perkins
 * @date   2013-11-14 11:35:49 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 26, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _vector_implementation_h_
#define _vector_implementation_h_


#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/utilities/complex.hpp>

namespace gridpack {
namespace math {

class ImplementationVisitor;
class ConstImplementationVisitor;

// -------------------------------------------------------------
//  class VectorImplementation
// -------------------------------------------------------------
class VectorImplementation 
  : private utility::Uncopyable,
    public parallel::Distributed
{
public:

  /// Default constructor.
  VectorImplementation(const parallel::Communicator& comm);

  /// Destructor
  ~VectorImplementation(void);

  /// Get the global length
  int size(void) const
  {
    return this->p_size();
  }

  /// Get the local length
  int localSize(void) const
  {
    return this->p_localSize();
  }

  /// Get the local min/max global indexes
  void localIndexRange(int& lo, int& hi) const
  {
    return this->p_localIndexRange(lo, hi);
  }

  /// Set an individual element
  void setElement(const int& i, const ComplexType& x)
  {
    this->p_setElement(i, x);
  }

  /// Set a range of elements (lo to hi-1)
  void setElementRange(const int& lo, const int& hi, ComplexType *x)
  {
    this->p_setElementRange(lo, hi, x);
  }

  /// Set an several elements
  void setElements(const int& n, const int *i, const ComplexType *x)
  {
    this->p_setElements(n, i, x);
  }

  /// Add to an individual element
  void addElement(const int& i, const ComplexType& x)
  {
    this->p_addElement(i, x);
  }

  /// Add to an several elements
  void addElements(const int& n, const int *i, const ComplexType *x)
  {
    this->p_addElements(n, i, x);
  }

  /// Get an individual element
  void getElement(const int& i, ComplexType& x) const
  {
    this->p_getElement(i, x);
  }

  /// Get an several elements
  void getElements(const int& n, const int *i, ComplexType *x) const
  {
    this->p_getElements(n, i, x);
  }

  /// Get a range of elements (lo to hi-1)
  void getElementRange(const int& lo, const int& hi, ComplexType *x) const
  {
    this->p_getElementRange(lo, hi, x);
  }

  /// Get all of vector elements (on all processes)
  void getAllElements(ComplexType *x) const
  {
    this->p_getAllElements(x);
  }


  /// Make all the elements zero
  void zero(void)
  {
    this->p_zero();
  }

  /// Make all the elements the specified value
  void fill(const ComplexType& v)
  {
    this->p_fill(v);
  }

  /// Compute the vector L1 norm (sum of absolute value)
  ComplexType norm1(void) const
  {
    return this->p_norm1();
  }

  /// Compute the vector L2 norm (root of sum of squares)
  ComplexType norm2(void) const
  {
    return this->p_norm2();
  }

  /// Compute the infinity (or maximum) norm
  ComplexType normInfinity(void) const
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

  // FIXME: more ...

  /// Make this instance ready to use
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

  // -------------------------------------------------------------
  // In-place Vector Operation Methods (change this instance)
  // -------------------------------------------------------------

  /// Make an exact replica of this instance
  VectorImplementation *clone(void) const
  {
    return this->p_clone();
  }

protected:

  /// Get the global vector length (specialized)
  virtual int p_size(void) const = 0;

  /// Get the size of the vector local part (specialized)
  virtual int p_localSize(void) const = 0;

  /// Get the local min/max global indexes (specialized)
  virtual void p_localIndexRange(int& lo, int& hi) const = 0;

  /// Set an individual element (specialized)
  virtual void p_setElement(const int& i, const ComplexType& x) = 0;

  /// Set an several elements (specialized)
  virtual void p_setElements(const int& n, const int *i, const ComplexType *x) = 0;

  /// Get a range of elements (lo to hi-1) (specialized)
  virtual void p_setElementRange(const int& lo, const int& hi, ComplexType *x);

  /// Add to an individual element (specialized)
  virtual void p_addElement(const int& i, const ComplexType& x) = 0;

  /// Add to an several elements (specialized)
  virtual void p_addElements(const int& n, const int *i, const ComplexType *x) = 0;

  /// Get an individual element (specialized)
  virtual void p_getElement(const int& i, ComplexType& x) const = 0;

  /// Get an several elements (specialized)
  virtual void p_getElements(const int& n, const int *i, ComplexType *x) const = 0;

  /// Get a range of elements (lo to hi-1) (specialized)
  virtual void p_getElementRange(const int& lo, const int& hi, ComplexType *x) const;

  /// Get all of vector elements (on all processes)
  virtual void p_getAllElements(ComplexType *x) const = 0;

  /// Make all the elements zero (specialized)
  virtual void p_zero(void) = 0;

  /// Fill all the elements with the specified value (specialized)
  virtual void p_fill(const ComplexType& v) = 0;

  /// Compute the vector L1 norm (sum of absolute value) (specialized)
  virtual ComplexType p_norm1(void) const = 0;

  /// Compute the vector L2 norm (root of sum of squares) (specialized)
  virtual ComplexType p_norm2(void) const = 0;

  /// Compute the vector infinity (or maximum) norm (specialized)
  virtual ComplexType p_normInfinity(void) const = 0;

  /// Replace all elements with its absolute value (specialized) 
  virtual void p_abs(void) = 0;

  /// Replace all elements with their real part (specialized)
  virtual void p_real(void);

  /// Replace all elements with their imaginary part (specialized)
  virtual void p_imaginary(void);

  /// Replace all elements with their complex conjugate
  virtual void p_conjugate(void);

  /// Make this instance ready to use
  virtual void p_ready(void) = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ImplementationVisitor& visitor) = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ConstImplementationVisitor& visitor) const = 0;

  /// Make an exact replica of this instance (specialized)
  virtual VectorImplementation *p_clone(void) const = 0;

  // -------------------------------------------------------------
  // In-place Vector Operation Methods (change this instance)
  // -------------------------------------------------------------

};

} // namespace math
} // namespace gridpack



#endif
