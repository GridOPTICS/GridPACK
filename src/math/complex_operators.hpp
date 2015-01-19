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
 * @file   complex_operators.hpp
 * @author William A. Perkins
 * @date   2015-01-19 10:52:09 d3g096
 * 
 * @brief This header provides type interregation utilities and some math
 * operators for the math library
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created January 16, 2015 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _complex_operators_hpp_
#define _complex_operators_hpp_

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include "gridpack/utilities/complex.hpp"

namespace gridpack {
namespace math {

/// Equate 
template <typename From, typename To>
inline void equate(To& t, const From& f)
{
  BOOST_STATIC_ASSERT(boost::is_same<From, To>::value);
  t = f;
}

template <>
inline void equate<RealType, ComplexType>(ComplexType& t, const RealType& f)
{
  t = f;
}

template <>
inline void equate<ComplexType, RealType>(RealType& t, const ComplexType& f)
{
  t = std::real(f);
}


/// Get the real part of a number
/** 
 * instantiate only for gridpack::RealType and gridpack::ComplexType;
 * other types should produce a compilation error (not a link error)
 * 
 * @param value 
 * 
 * @return 
 */
template <typename T> inline RealType realpart(const T& value);


// Instantiation for real values
template <> 
inline RealType realpart<RealType>(const RealType& value)
{
  return value;
}

// Instantiation for complex values
template <> 
inline RealType realpart<ComplexType>(const ComplexType& value)
{
  return std::real(value);
}

  

} // namespace math
} // namespace gridpack
#endif
