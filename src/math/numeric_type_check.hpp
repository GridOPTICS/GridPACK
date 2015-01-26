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
 * @file   numeric_type_check.hpp
 * @author William A. Perkins
 * @date   2015-01-22 15:06:17 d3g096
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


#ifndef _numeric_type_check_hpp_
#define _numeric_type_check_hpp_

#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>

#include <gridpack/utilities/complex.hpp>

namespace gridpack {
namespace math {

/// This is used to make sure the numeric type is supported
/**
 * GridPACK only supports a limited number of element types for
 * matrices and vectors, specifically @c double and @c
 * "std::complex<double>. This will cause a compilation failure if the
 * type is not supported.
 * 
 */
template <typename T>
struct TypeCheck {
private:
  typedef boost::mpl::vector<RealType, ComplexType> Allowable;
  typedef boost::mpl::end<Allowable>::type NotFound;
  typedef typename boost::mpl::find<Allowable, T>::type Found;
public:
  typedef boost::mpl::not_<
    boost::is_same<Found, NotFound>
  > OK;
  static const bool check = OK::value;
  static const bool isComplex = boost::is_same<T, ComplexType>::value;
  template <typename Other>
  struct isSame : boost::is_same<T, Other>::type
  {};
};

} // namespace math
} // namespace gridpack




#endif
