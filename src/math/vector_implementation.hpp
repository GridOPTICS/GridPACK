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
 * @date   2014-10-21 09:23:22 d3g096
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
  VectorImplementation(const parallel::Communicator& comm);

  /// Destructor
  ~VectorImplementation(void);

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

} // namespace math
} // namespace gridpack



#endif
