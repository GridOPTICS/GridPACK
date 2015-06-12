// -------------------------------------------------------------
// file: vector_implementation.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 14, 2014 by William A. Perkins
// Last Change: 2015-06-12 09:11:48 d3g096
// -------------------------------------------------------------


#include "vector_implementation.hpp"

namespace gridpack {
namespace math {

template <>
void
VectorImplementation<ComplexType>::p_real(void)
{
  int lo, hi;
  this->localIndexRange(lo, hi);
  std::vector<TheType> x(hi-lo);
  this->getElementRange(lo, hi, &x[0]);
  std::vector<TheType>::iterator i;
  for (i = x.begin(); i != x.end(); ++i) {
    *i = std::real(*i);
  }
  this->setElementRange(lo, hi, &x[0]);
  this->ready();
}

template <>
void
VectorImplementation<RealType>::p_real(void)
{
  // do nothing
}

// -------------------------------------------------------------
// VectorImplementation::p_imaginary
// -------------------------------------------------------------
template <>
void
VectorImplementation<ComplexType>::p_imaginary(void)
{
  int lo, hi;
  this->localIndexRange(lo, hi);
  std::vector<TheType> x(hi-lo);
  this->getElementRange(lo, hi, &x[0]);
  for (std::vector<TheType>::iterator i = x.begin();
       i != x.end(); ++i) {
    *i = imag(*i);
  }
  this->setElementRange(lo, hi, &x[0]);
  this->ready();
}

template <>
void
VectorImplementation<RealType>::p_imaginary(void)
{
  this->p_zero();
}

} // namespace math
} // namespace gridpack

