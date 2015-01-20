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
 * @file   value_transfer.hpp
 * @author William A. Perkins
 * @date   2015-01-20 10:55:27 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created January 20, 2015 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _value_transfer_hpp_
#define _value_transfer_hpp_

#include <boost/shared_array.hpp>
#include <boost/static_assert.hpp>

#include "numeric_type_check.hpp"

#include "gridpack/utilities/uncopyable.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class valueTransfer
// -------------------------------------------------------------
/**
 * This provides a consistent way to 
 * 
 */

template <typename FromType, typename ToType>
class ValueTransfer 
  : private utility::Uncopyable
{
private:

  BOOST_STATIC_ASSERT(TypeCheck<FromType>::OK::value);
  BOOST_STATIC_ASSERT(TypeCheck<ToType>::OK::value);

  /// A thing to provide a null delete operation
  struct null_deleter {
    void operator()(void const *p) const {}
  };
  
  
public:

  /// Default constructor.
  ValueTransfer(const unsigned int& from_size, FromType* from)
    : utility::Uncopyable(),
      p_fromSize(from_size), p_from(from)
  {
      p_setup();
  }

  /// Destructor
  ~ValueTransfer(void) 
  {  }

  /// Get the length of the "to" buffer
  unsigned int size(void) const
  {
    return p_toSize;
  }

  /// Get a pointer to the "to" buffer
  ToType *to(void)
  {
    return p_to.get();
  }

protected:

  /// The size of the from buffer
  const unsigned int p_fromSize;

  /// The from buffer
  FromType *p_from;

  /// The size of the to buffer
  unsigned int p_toSize;

  /// The "to" buffer
  boost::shared_array<ToType> p_to;

  /// Do the setup  
  void p_setup(void);

};

template <>
void
ValueTransfer<RealType, RealType>::p_setup(void)
{
  p_toSize = p_fromSize;
  p_to.reset(p_from, null_deleter());
}

template <>
void
ValueTransfer<ComplexType, ComplexType>::p_setup(void)
{
  p_toSize = p_fromSize;
  p_to.reset(p_from, null_deleter());
}

template <>
void
ValueTransfer<RealType, ComplexType>::p_setup(void)
{
  const unsigned int TWO(2);
  BOOST_ASSERT(p_fromSize > 1);
  p_toSize = p_fromSize/TWO;
  p_to.reset(new ComplexType[p_toSize]);
  unsigned int fidx(0), tidx(0);
  for (fidx = 0; fidx < p_fromSize; fidx += TWO, ++tidx) {
    p_to.get()[tidx] = ComplexType(p_from[fidx], p_from[fidx+1]);
  }
}

template<>
void
ValueTransfer<ComplexType, RealType>::p_setup(void)
{
  const unsigned int TWO(2);
  p_toSize = p_fromSize*TWO;
  p_to.reset(new RealType[p_toSize]);
  unsigned int fidx(0), tidx(0);
  for (fidx = 0; fidx < p_fromSize; ++fidx, ++tidx += TWO) {
    p_to.get()[tidx] = std::real(p_from[fidx]);
    p_to.get()[tidx+1] = std::imag(p_from[fidx]);
  }
}
  
} // namespace math
} // namespace gridpack

#endif
