// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   vector_implementation.cpp
 * @author William A. Perkins
 * @date   2014-10-30 09:29:47 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/utilities/complex.hpp"

inline 
gridpack::ComplexType
stupid_exp(const gridpack::ComplexType& x)
{
  return std::exp(x);
}



#include "vector_implementation.hpp"
#include <cmath>
#include <vector>
#include <boost/iterator/counting_iterator.hpp>

namespace gridpack {
namespace math {


// -------------------------------------------------------------
//  class VectorImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// VectorImplementation::p_exp
// -------------------------------------------------------------
inline void
VectorImplementation::p_exp(void)
{
  int lo, hi;
  this->localIndexRange(lo, hi);
  std::vector<ComplexType> x(hi-lo);
  this->getElementRange(lo, hi, &x[0]);
  for (std::vector<ComplexType>::iterator i = x.begin();
       i != x.end(); ++i) {
    ComplexType x(*i);
    *i = stupid_exp(x);
  }
  this->setElementRange(lo, hi, &x[0]);
  this->ready();
}

// -------------------------------------------------------------
// VectorImplementation::p_conjugate
// -------------------------------------------------------------
inline void
VectorImplementation::p_conjugate(void)
{
  int lo, hi;
  this->localIndexRange(lo, hi);
  std::vector<ComplexType> x(hi-lo);
  this->getElementRange(lo, hi, &x[0]);
  for (std::vector<ComplexType>::iterator i = x.begin();
       i != x.end(); ++i) {
    *i = conj(*i);
  }
  this->setElementRange(lo, hi, &x[0]);
  this->ready();
}


} // namespace math
} // namespace gridpack

