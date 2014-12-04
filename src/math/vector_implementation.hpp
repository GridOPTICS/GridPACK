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
 * @date   2014-12-04 10:46:00 d3g096
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
#include <gridpack/math/vector_interface.hpp>

namespace gridpack {
namespace math {

class ImplementationVisitor;
class ConstImplementationVisitor;

// -------------------------------------------------------------
//  class VectorImplementation
// -------------------------------------------------------------
class VectorImplementation 
  : private utility::Uncopyable,
    public parallel::Distributed,
    public BaseVectorInterface<ComplexType>
{
public:

  /// Default constructor.
  VectorImplementation(const parallel::Communicator& comm)
    : utility::Uncopyable(), parallel::Distributed(comm)
  {}

  /// Destructor
  ~VectorImplementation(void)
  {}

  // -------------------------------------------------------------
  // Allow visitors
  // -------------------------------------------------------------

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

  /// Get a range of elements (lo to hi-1) (specialized)
  void p_setElementRange(const IdxType& lo, const IdxType& hi, TheType *x);

  /// Get a range of elements (lo to hi-1) (specialized)
  void p_getElementRange(const IdxType& lo, const IdxType& hi, TheType *x) const;

  /// Replace all elements with their real part (specialized)
  void p_real(void);

  /// Replace all elements with their imaginary part (specialized)
  void p_imaginary(void);

  /// Replace all elements with their complex conjugate
  void p_conjugate(void);

  /// Replace all elements with its exponential (specialized)
  void p_exp(void);

  /// Allow visits by implementation visitors
  virtual void p_accept(ImplementationVisitor& visitor) = 0;
  
  /// Allow visits by implementation visitors
  virtual void p_accept(ConstImplementationVisitor& visitor) const = 0;

  /// Make an exact replica of this instance (specialized)
  virtual VectorImplementation *p_clone(void) const = 0;
  
};

// -------------------------------------------------------------
// VectorImplementation::p_setElementRange
// -------------------------------------------------------------
inline void
VectorImplementation::p_setElementRange(const int& lo, const int& hi, ComplexType *x)
{
  std::vector<int> i;
  i.reserve(hi-lo);
  std::copy(boost::counting_iterator<int>(lo),
            boost::counting_iterator<int>(hi),
            std::back_inserter(i));
  this->p_setElements(i.size(), &i[0], x);
}

// -------------------------------------------------------------
// VectorImplementation::p_getElementRange
// -------------------------------------------------------------
inline void
VectorImplementation::p_getElementRange(const int& lo, const int& hi, ComplexType *x) const
{
  std::vector<int> i;
  i.reserve(hi-lo);
  std::copy(boost::counting_iterator<int>(lo),
            boost::counting_iterator<int>(hi),
            std::back_inserter(i));
  this->p_getElements(i.size(), &i[0], x);
}

// -------------------------------------------------------------
// VectorImplementation::p_real
// -------------------------------------------------------------
inline void
VectorImplementation::p_real(void)
{
  int lo, hi;
  this->localIndexRange(lo, hi);
  std::vector<ComplexType> x(hi-lo);
  this->getElementRange(lo, hi, &x[0]);
  for (std::vector<ComplexType>::iterator i = x.begin();
       i != x.end(); ++i) {
    *i = std::real(*i);
  }
  this->setElementRange(lo, hi, &x[0]);
  this->ready();
}
  
// -------------------------------------------------------------
// VectorImplementation::p_imaginary
// -------------------------------------------------------------
inline void
VectorImplementation::p_imaginary(void)
{
  int lo, hi;
  this->localIndexRange(lo, hi);
  std::vector<ComplexType> x(hi-lo);
  this->getElementRange(lo, hi, &x[0]);
  for (std::vector<ComplexType>::iterator i = x.begin();
       i != x.end(); ++i) {
    *i = imag(*i);
  }
  this->setElementRange(lo, hi, &x[0]);
  this->ready();
}

} // namespace math
} // namespace gridpack



#endif
