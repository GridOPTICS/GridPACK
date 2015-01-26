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
 * @date   2015-01-22 15:25:09 d3g096
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
 * This provides a consistent way to transfer values between RealType
 * and ComplexType arrays, being (hopefully) smart about when memory
 * is allocated. The transfer works both ways: complex <-> real. This
 * is not approprate for the case of real math with a complex-based
 * math library!
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
  ValueTransfer(const unsigned int& from_size, FromType* from, ToType* to = NULL)
    : utility::Uncopyable(),
      p_fromSize(from_size), p_from(from),
      p_toSize(), p_to()
  {
    p_toSize = p_computeToSize();
    BOOST_ASSERT(p_toSize > 0);
    p_setup(to);
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

  /// Two
  static const unsigned int TWO = 2;

  /// The size of the from buffer
  const unsigned int p_fromSize;

  /// The from buffer
  FromType *p_from;

  /// The size of the to buffer
  unsigned int p_toSize;

  /// The "to" buffer
  boost::shared_array<ToType> p_to;

  /// Compute the size of the "to" buffer
  virtual inline unsigned int p_computeToSize(void);

  /// Do the setup  
  inline void p_setup(ToType *to);

  /// Copy the value
  virtual inline void p_copy(void);

};


// -------------------------------------------------------------
// ValueTransfer<>::p_computeToSize
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline unsigned int
ValueTransfer<FromType, ToType>::p_computeToSize(void)
{
  BOOST_STATIC_ASSERT(boost::is_same<FromType, ToType>::value);
  return p_fromSize;
}

template <>
inline unsigned int
ValueTransfer<RealType, ComplexType>::p_computeToSize(void)
{
  return p_fromSize/TWO;
};

template <>
inline unsigned int
ValueTransfer<ComplexType, RealType>::p_computeToSize(void)
{
  return p_fromSize*TWO;
};

// -------------------------------------------------------------
// ValueTransfer<>::p_setup
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline void  
ValueTransfer<FromType, ToType>::p_setup(ToType *to)
{
  if (to != NULL) {
    p_to.reset(to, null_deleter());
  } else {
    p_to.reset(new ToType[p_toSize]);
  }
  p_copy();
}

template <>
inline void  
ValueTransfer<RealType, RealType>::p_setup(RealType *to)
{
  if (to != NULL) {
    p_to.reset(to, null_deleter());
    p_copy();
  } else {
    p_to.reset(p_from, null_deleter());
  }
};

template <>
inline void  
ValueTransfer<ComplexType, ComplexType>::p_setup(ComplexType *to)
{
  if (to != NULL) {
    p_to.reset(to, null_deleter());
    p_copy();
  } else {
    p_to.reset(p_from, null_deleter());
  }
};



// -------------------------------------------------------------
// ValueTransfer<>::p_copy
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline void
ValueTransfer<FromType, ToType>::p_copy(void)
{
  BOOST_STATIC_ASSERT(boost::is_same<FromType, ToType>::value);
  if (p_from != p_to.get()) 
    std::copy(p_from, p_from + p_fromSize, p_to.get());
}

template <>
inline void
ValueTransfer<RealType, ComplexType>::p_copy(void)
{
  unsigned int fidx(0), tidx(0);
  for (fidx = 0; fidx < p_fromSize; fidx += TWO, ++tidx) {
    p_to.get()[tidx] = ComplexType(p_from[fidx], p_from[fidx+1]);
  }
};

template <>
inline void
ValueTransfer<ComplexType, RealType>::p_copy(void)
{
  unsigned int fidx(0), tidx(0);
  for (fidx = 0; fidx < p_fromSize; ++fidx, ++tidx += TWO) {
    p_to.get()[tidx] = std::real(p_from[fidx]);
    p_to.get()[tidx+1] = std::imag(p_from[fidx]);
  }
}

// -------------------------------------------------------------
// storage_size
// -------------------------------------------------------------
/// Get the number of library elements used to represent a scalar
/** 
 * complex -> complex -> 1
 * real -> real -> 1
 * real -> complex -> 1
 * complex -> real -> 2
 * 
 * @return 
 */
template <typename ScalarType, typename LibraryType>
struct storage_size
  : boost::mpl::if_<
      typename boost::mpl::and_< 
        boost::is_same<ScalarType, ComplexType>,
        boost::is_same<LibraryType, RealType> >::type,
      boost::mpl::int_<2>,
    boost::mpl::int_<1> 
  >::type
{};



} // namespace math
} // namespace gridpack



#endif
