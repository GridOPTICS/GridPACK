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
 * @date   2015-06-12 09:12:12 d3g096
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
/// 
/** 
 * This allows a general way to do @c "x = y" 
 * 
 * @param f 
 * 
 * @return 
 */
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
/** 
 * This function applies an operator to each element in an array of
 * type T that is stored in an array of type @c StorageType. The
 * values array @c are replaced.
 *
 * In this implementation, it is assumed that types @c T and @c
 * StorageType are the same or castable to each other.
 * 
 * @param size the number of type @c T elements to be operated on
 * @param x the elements stored as @c StorageType
 * @param op 
 */
template <typename T, typename StorageType>
inline void
unary_operation(const unsigned int& size, StorageType *x, 
                base_unary_function<T>& op)
{
  for (unsigned int i = 0; i < size; ++i) {
    x[i] = op(x[i]);
  }
}


/** 
 * In this specialization, ComplexType values are stored in a
 * RealType vector, two RealType values are used for each ComplexType.
 * So, a temporary ComplexType is made from 2 RealType entries,
 * operated on, then split back into two RealType elements.
 * 
 * @param size 
 * @param x 
 * @param op 
 * 
 * @return 
 */
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

/** 
 * In this specialization, RealType values are stored in a ComplexType
 * array. It is assumed that only the real part is used.
 * 
 * @param size 
 * @param x 
 * @param op 
 * 
 * @return 
 */
template <>
inline void
unary_operation<RealType, ComplexType>(const unsigned int& size, ComplexType *x, 
                                       base_unary_function<RealType>& op)
{
  for (unsigned int i = 0; i < size; ++i) {
    RealType tmp(std::real(x[i]));
    x[i] = op(tmp);
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
// addvalue
// -------------------------------------------------------------
template <typename T> 
struct addvalue : public base_unary_function<T>
{
  const T value;
  addvalue(const T& v) : value(v) {}
  T operator() (const T& x) const { return x+value; }
};


// -------------------------------------------------------------
// multiply
// -------------------------------------------------------------
template <typename T> 
struct multiplyvalue : public base_unary_function<T>
{
  const T factor;
  multiplyvalue(const T& f) : factor(f) {}
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
    return one/x; 
  }
};

// -------------------------------------------------------------
// conjugate
// -------------------------------------------------------------
template <typename T> 
struct conjugate_value : public base_unary_function<T>
{
  inline T operator() (const T& x) const;
};

template <>
inline RealType
conjugate_value<RealType>::operator() (const RealType& x) const
{
  return x;
}

template <>
inline ComplexType
conjugate_value<ComplexType>::operator() (const ComplexType& x) const
{
  return std::conj(x);
}

// -------------------------------------------------------------
// real_value
// -------------------------------------------------------------
template <typename T> 
struct real_value : public base_unary_function<T>
{
  inline T operator() (const T& x) const;
};

template <>
inline RealType
real_value<RealType>::operator() (const RealType& x) const
{
  return x;
}

template <>
inline ComplexType
real_value<ComplexType>::operator() (const ComplexType& x) const
{
  return std::real(x);
}

// -------------------------------------------------------------
// imaginary_value
// -------------------------------------------------------------
template <typename T> 
struct imaginary_value : public base_unary_function<T>
{
  inline T operator() (const T& x) const;
};

template <>
inline RealType
imaginary_value<RealType>::operator() (const RealType& x) const
{
  return 0.0;
}

template <>
inline ComplexType
imaginary_value<ComplexType>::operator() (const ComplexType& x) const
{
  return std::imag(x);
}


// -------------------------------------------------------------
// base_accumulator_function
// -------------------------------------------------------------
template <typename T, typename ResultType>
struct base_accumulator_function : public std::unary_function<T, void>
{
  ResultType accum;
  virtual void operator() (const T& t) = 0;
  virtual ResultType result(void) const 
  {
    return accum;
  }
  base_accumulator_function(void) 
    : std::unary_function<T, void>(), accum(0.0)
  {}
};

// -------------------------------------------------------------
// accumulator_operation
// -------------------------------------------------------------

template <typename T, typename StorageType>
inline void
accumulator_operation(const unsigned int& size, const StorageType *x, 
                      base_accumulator_function<T, RealType>& op)
{
  for (unsigned int i = 0; i < size; ++i) {
    op(x[i]);
  }
}

template <>
inline void
accumulator_operation<RealType, ComplexType>(const unsigned int& size, 
                                             const ComplexType *x, 
                                             base_accumulator_function<RealType, RealType>& op)
{
  for (unsigned int i = 0; i < size; i += 1) {
    RealType tmp(real(x[i]));
    op(tmp);
  }
}

template <>
inline void
accumulator_operation<ComplexType, RealType>(const unsigned int& size, 
                                             const RealType *x, 
                                             base_accumulator_function<ComplexType, RealType>& op)
{
  for (unsigned int i = 0; i < size; i += 2) {
    ComplexType tmp(x[i], x[i+1]);
    op(tmp);
  }
}
  

// -------------------------------------------------------------
// l1_norm
// -------------------------------------------------------------
template <typename T>
struct l1_norm : public base_accumulator_function<T, RealType>
{
  inline void operator() (const T& t);
};

template <>
inline void
l1_norm<RealType>::operator() (const RealType& x)
{
  accum += ::fabs(x);
}

template <>
inline void
l1_norm<ComplexType>::operator() (const ComplexType& x)
{
  accum += std::abs(x);
}

// -------------------------------------------------------------
// l2_norm
// -------------------------------------------------------------
template <typename T>
struct l2_norm : public base_accumulator_function<T, RealType>
{
  inline void operator() (const T& x);
};

template <>
inline void
l2_norm<RealType>::operator() (const RealType& x)
{
  accum += x*x;
}

template <>
inline void
l2_norm<ComplexType>::operator() (const ComplexType& x)
{
  RealType rx(std::abs(x));
  accum += rx*rx;
}

// -------------------------------------------------------------
// infinity_norm
// -------------------------------------------------------------
template <typename T>
struct infinity_norm : public base_accumulator_function<T, RealType>
{
  inline void operator() (const T& x);
};

template <>
inline void
infinity_norm<RealType>::operator() (const RealType& x)
{
  
  accum = std::max(accum, ::fabs(x));
}

template <>
inline void
infinity_norm<ComplexType>::operator() (const ComplexType& x)
{
  
  accum = std::max(accum, std::abs(x));
}


// -------------------------------------------------------------
// binary_operation
// -------------------------------------------------------------
template <typename T> 
struct base_binary_function : public std::binary_function<T, T, T>
{
  virtual T operator() (const T& x1, const T& x2) const = 0;
};

// -------------------------------------------------------------
// binary_operation
// -------------------------------------------------------------
/** 
 * This function combines two arrays of type @c T that are stored as
 * type @c StorageType using some operation.  The values in @c x1 are
 * replaced with the result of the operation.
 *
 * In this implementation, it is assumed that types @c T and @c
 * StorageType are the same or castable to each other and that one @c
 * StorageType element is used to represent one @c T element.
 * 
 * @param size 
 * @param x1 
 * @param x2 
 * @param op 
 */
template <typename T, typename StorageType>
inline void
binary_operation(const unsigned int& size, 
                 StorageType *x1, const StorageType *x2,
                 base_binary_function<T>& op)
{
  for (unsigned int i = 0; i < size; i += 1) {
    T result = op(x1[i], x2[i]);
    x1[i] = result;
  }
}


/** 
 * This function combines two arrays of type @c T that are stored as
 * type @c StorageType using some operation.  The values in @c x1 are
 * replaced with the result of the operation.
 *
 * In this specialization, ComplexType values are stored in a
 * RealType vector, two RealType values are used for each ComplexType.
 * So, temporary ComplexType's are made from 2 RealType entries,
 * operated on, then split back into two RealType elements.
 * 
 * @param size 
 * @param x1 
 * @param x2 
 * @param op 
 * 
 * @return 
 */
template <>
inline void
binary_operation<ComplexType, RealType>(const unsigned int& size, 
                                        RealType *x1, const RealType *x2,
                                        base_binary_function<ComplexType>& op)
{
  for (unsigned int i = 0; i < size; i += 2) {
    ComplexType tmp1(x1[i], x1[i+1]);
    ComplexType tmp2(x2[i], x2[i+1]);
    ComplexType result(op(tmp1, tmp2));
    x1[i] = std::real(result);
    x1[i+1] = std::imag(result);
  }
}

/** 
 * This function combines two arrays of type @c T that are stored as
 * type @c StorageType using some operation.  The values in @c x1 are
 * replaced with the result of the operation.
 *
 * 
 * @param size 
 * @param x1 
 * @param x2 
 * @param op 
 * 
 * @return 
 */
template <>
inline void
binary_operation<RealType, ComplexType>(const unsigned int& size, 
                                        ComplexType *x1, const ComplexType *x2,
                                        base_binary_function<RealType>& op)
{
  for (unsigned int i = 0; i < size; ++i) {
    RealType tmp1(std::real(x1[i]));
    RealType tmp2(std::real(x2[i]));
    RealType result(op(tmp1, tmp2));
    x1[i] = result;
  }
}

// -------------------------------------------------------------
// multiplyvalue2
// -------------------------------------------------------------
template <typename T> 
struct multiplyvalue2 : public base_binary_function<T>
{
  T operator() (const T& x1, const T& x2) const { return x1*x2; }
};

// -------------------------------------------------------------
// dividevalue2
// -------------------------------------------------------------
template <typename T> 
struct dividevalue2 : public base_binary_function<T>
{
  T operator() (const T& x1, const T& x2) const { return x1/x2; }
};

// -------------------------------------------------------------
// scaleAdd2
// -------------------------------------------------------------
template <typename T>
struct scaleAdd2 : public base_binary_function<T>
{
  T alpha;
  T operator() (const T& x1, const T& x2) const { return x1 + alpha*x2; }
  scaleAdd2(const T& a = 1.0)
    : base_binary_function<T>(), alpha(a)
  {}
};


} // namespace math
} // namespace gridpack
#endif
