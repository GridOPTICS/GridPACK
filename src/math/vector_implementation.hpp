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
 * @date   2019-11-20 10:01:00 d3g096
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
#include <gridpack/utilities/exception.hpp>
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

  /// A place to store indexes for set/get element range
  mutable std::vector<IdxType> p_rangeIdx;

  /// A place to store local elements, if necessary
  boost::scoped_array<TheType> p_localElements;

  /// Make ::p_rangeIdx
  void p_buildRangeIdx(void) const
  {
    if (p_rangeIdx.empty()) {
      IdxType n(this->size());
      p_rangeIdx.reserve(n);
      std::copy(boost::counting_iterator<int>(0),
                boost::counting_iterator<int>(n),
                std::back_inserter(p_rangeIdx));
    }
  }

  /// 
  virtual void p_applyOperation(base_unary_function<TheType>& op) = 0;


  /// Set a range of elements (lo to hi-1) (specialized)
  void p_setElementRange(const IdxType& lo, const IdxType& hi, TheType *x)
  {
    p_buildRangeIdx();
    IdxType n(hi - lo);
    this->p_setElements(n, &p_rangeIdx[lo], x);
  }

  /// Set a range of elements (lo to hi-1) (specialized)
  void p_addElementRange(const IdxType& lo, const IdxType& hi, TheType *x)
  {
    p_buildRangeIdx();
    IdxType n(hi - lo);
    this->p_addElements(n, &p_rangeIdx[lo], x);
  }

  /// Get a range of elements (lo to hi-1) (specialized)
  void p_getElementRange(const IdxType& lo, const IdxType& hi, TheType *x) const
  {
    p_buildRangeIdx();
    IdxType n(hi - lo);
    this->p_getElements(n, &p_rangeIdx[lo], x);
  }

  /// Get an array of local elements (specialized)
  /** 
   * A very generic implementation. 
   * 
   * 
   * @return 
   */
  TheType *p_getLocalElements(void)
  {
    IdxType n(this->localSize());
    IdxType lo, hi;
    this->localIndexRange(lo, hi);

    if (!p_localElements) {
      p_localElements.reset(new TheType[n]);
    }
    this->getElementRange(lo, hi, p_localElements.get());
    return p_localElements.get();
  }

  /// release the local elements array produced by getLocalElements() (specialized)
  void p_releaseLocalElements(TheType *array)
  {
    // some things that should not happen
    
    if (!p_localElements) {
      throw gridpack::Exception("VectorImplementation: "
                                "releaseLocalElements called w/o getLocalElements");
    }
    if (array != p_localElements.get()) {
      throw gridpack::Exception("VectorImplementation: "
                                "releaseLocalElements array does not match getLocalElements");
    }

    IdxType lo, hi;
    this->localIndexRange(lo, hi);
    this->setElementRange(lo, hi, p_localElements.get());
    // may want to keep this around if operation is done alot
    p_localElements.reset();
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
