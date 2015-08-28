// -------------------------------------------------------------
// file: expression.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 28, 2015 by William A. Perkins
// Last Change: 2015-08-28 16:40:39 d3g096
// -------------------------------------------------------------

#include <boost/assert.hpp>
#include <boost/format.hpp>
#include "expression.hpp"

namespace gridpack {
namespace optimization {

template <typename T>
std::string
ConstantExpression<T>::p_render(void) const
{
  BOOST_ASSERT(false);
}

template <>
std::string
ConstantExpression<int>::p_render(void) const
{
  return boost::str(boost::format("%d") % p_value);
}

template <>
std::string
ConstantExpression<double>::p_render(void) const
{
  return boost::str(boost::format("%g") % p_value);
}


} // namespace optimization
} // namespace gridpack

