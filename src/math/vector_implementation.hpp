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
 * @date   2015-01-26 10:50:14 d3g096
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
#include <gridpack/math/complex_operators.hpp>

namespace gridpack {
namespace math {

class ImplementationVisitor;
class ConstImplementationVisitor;

// -------------------------------------------------------------
//  class VectorImplementation
// -------------------------------------------------------------
template <typename T, typename I = int>
class VectorImplementation 
  : private utility::Uncopyable,
    public parallel::Distributed,
    public BaseVectorInterface<T, I>
{
public:

  typedef typename BaseVectorInterface<T, I>::IdxType IdxType;
  typedef typename BaseVectorInterface<T, I>::TheType TheType;

  /// Default constructor.
  VectorImplementation(const parallel::Communicator& comm)
    : utility::Uncopyable(), parallel::Distributed(comm)
  {}

  /// Destructor
  ~VectorImplementation(void)
  {}

  /// Apply some operator on the vector elements
  void applyOperation(base_unary_function<TheType>& op)
  {
    this->p_applyOperation(op);
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

  /// 
  virtual void p_applyOperation(base_unary_function<TheType>& op) = 0;


  /// Get a range of elements (lo to hi-1) (specialized)
  void p_setElementRange(const IdxType& lo, const IdxType& hi, TheType *x)
  {
    std::vector<int> i;
    i.reserve(hi-lo);
    std::copy(boost::counting_iterator<int>(lo),
              boost::counting_iterator<int>(hi),
              std::back_inserter(i));
    this->p_setElements(i.size(), &i[0], x);
  }

  /// Get a range of elements (lo to hi-1) (specialized)
  void p_getElementRange(const IdxType& lo, const IdxType& hi, TheType *x) const
  {
    std::vector<int> i;
    i.reserve(hi-lo);
    std::copy(boost::counting_iterator<int>(lo),
              boost::counting_iterator<int>(hi),
              std::back_inserter(i));
    this->p_getElements(i.size(), &i[0], x);
  }

  /// Replace all elements with their real part (specialized)
  void p_real(void);

  /// Replace all elements with their imaginary part (specialized)
  void p_imaginary(void);

  /// Make an exact replica of this instance (specialized)
  virtual VectorImplementation *p_clone(void) const = 0;
  
};


} // namespace math
} // namespace gridpack



#endif
