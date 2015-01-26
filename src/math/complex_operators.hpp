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
 * @date   2015-01-26 09:57:23 d3g096
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

#include <functional>
#include <algorithm>

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include "gridpack/utilities/complex.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
// equate
// -------------------------------------------------------------
template <typename To, typename From>
inline To 
equate(const From& f)
{
  To t;
  t = f;
  return t;
}

template <>
inline RealType
equate<RealType, ComplexType>(const ComplexType& f)
{
  RealType t;
  t = std::real(f);
  return t;
}
  
/*

template <typename To, typename From>
inline void equate(To& t, const From& f)
{
  BOOST_STATIC_ASSERT(boost::is_same<From, To>::value);
  t = f;
}

template <>
inline void equate<ComplexType, RealType>(ComplexType& t, const RealType& f)
{
  t = f;
}

template <>
inline void equate<RealType, ComplexType>(RealType& t, const ComplexType& f)
{
  t = std::real(f);
}
*/

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

// -------------------------------------------------------------
// base_unary_function
// -------------------------------------------------------------
template <typename T>
struct base_unary_function : public std::unary_function<T, T>
{
  virtual T operator() (const T& t) const = 0;
};

// -------------------------------------------------------------
// unary_operation
// -------------------------------------------------------------
template <typename T, typename StorageType>
inline void
unary_operation(const unsigned int& size, StorageType *x, 
                base_unary_function<T>& op)
{
  for (unsigned int i = 0; i < size; i += 2) {
    x[i] = op(x[i]);
  }
}

template <>
inline void
unary_operation<ComplexType, RealType>(const unsigned int& size, RealType *x, 
                                       base_unary_function<ComplexType>& op)
{
  for (unsigned int i = 0; i < size; i += 2) {
    ComplexType tmp(x[i], x[i+1]);
    tmp = op(tmp);
    x[i] = std::real(tmp);
    x[i+1] = std::imag(tmp);
  }
}

template <>
inline void
unary_operation<RealType, ComplexType>(const unsigned int& size, ComplexType *x, 
                                       base_unary_function<RealType>& op)
{
  for (unsigned int i = 0; i < size; ++i) {
    ComplexType tmp(x[i]);
    x[i] = op(std::real(tmp));
  }
}

// -------------------------------------------------------------
// setvalue
// -------------------------------------------------------------
template <typename T> 
struct setvalue : public base_unary_function<T>
{
  const T value;
  setvalue(const T& v) : value(v) {}
  T operator() (const T& x) const { return value; }
};


// -------------------------------------------------------------
// multiply
// -------------------------------------------------------------
template <typename T> 
struct multiply : public base_unary_function<T>
{
  const T factor;
  multiply(const T& f) : factor(f) {}
  T operator() (const T& x) const { return x*factor; }
};

// -------------------------------------------------------------
// absvalue
// -------------------------------------------------------------
template <typename T> 
struct absvalue : public base_unary_function<T>
{
  T operator() (const T& x) const { 
    return x; 
  }
};

template <> 
struct absvalue<RealType> : public base_unary_function<RealType>
{
  RealType operator() (const RealType& x) const { 
    return fabs(x);
  }
};

template <> 
struct absvalue<ComplexType> : public base_unary_function<ComplexType>
{
  ComplexType operator() (const ComplexType& x) const { 
    return std::abs(x);
  }
};


// -------------------------------------------------------------
// exponential
// -------------------------------------------------------------
template <typename T> 
struct exponential : public base_unary_function<T>
{
  T operator() (const T& x) const { return exp(x); }
};

// -------------------------------------------------------------
// reciprocal
// -------------------------------------------------------------
template <typename T> 
struct reciprocal : public base_unary_function<T>
{
  T operator() (const T& x) const { 
    T one(1.0);
    return 1.0/x; 
  }
};


  

} // namespace math
} // namespace gridpack
#endif
