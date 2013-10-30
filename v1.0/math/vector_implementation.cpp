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
 * @date   2013-10-28 13:17:25 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <boost/iterator/counting_iterator.hpp>

#include "vector_implementation.hpp"

namespace gridpack {
namespace math {


// -------------------------------------------------------------
//  class VectorImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// VectorImplementation:: constructors / destructor
// -------------------------------------------------------------
VectorImplementation::VectorImplementation(const parallel::Communicator& comm)
  : utility::Uncopyable(), parallel::Distributed(comm)
{
  
}

VectorImplementation::~VectorImplementation(void)
{
}

// -------------------------------------------------------------
// VectorImplementation::p_setElementRange
// -------------------------------------------------------------
void
VectorImplementation::p_setElementRange(const int& lo, const int& hi, ComplexType *x)
{
  std::vector<int> i;
  i.reserve(hi-lo);
  std::copy(boost::counting_iterator<int>(lo),
            boost::counting_iterator<int>(hi),
            std::back_inserter(i));
  this->p_setElements(i.size(), &i[0], x);
}

// -------------------------------------------------------------
// VectorImplementation::p_getElementRange
// -------------------------------------------------------------
void
VectorImplementation::p_getElementRange(const int& lo, const int& hi, ComplexType *x) const
{
  std::vector<int> i;
  i.reserve(hi-lo);
  std::copy(boost::counting_iterator<int>(lo),
            boost::counting_iterator<int>(hi),
            std::back_inserter(i));
  this->p_getElements(i.size(), &i[0], x);
}

// -------------------------------------------------------------
// VectorImplementation::p_real
// -------------------------------------------------------------
void
VectorImplementation::p_real(void)
{
  int lo, hi;
  this->localIndexRange(lo, hi);
  std::vector<ComplexType> x(hi-lo);
  this->getElementRange(lo, hi, &x[0]);
  for (std::vector<ComplexType>::iterator i = x.begin();
       i != x.end(); ++i) {
    *i = std::real(*i);
  }
  this->setElementRange(lo, hi, &x[0]);
  this->ready();
}
  
// -------------------------------------------------------------
// VectorImplementation::p_imaginary
// -------------------------------------------------------------
void
VectorImplementation::p_imaginary(void)
{
  int lo, hi;
  this->localIndexRange(lo, hi);
  std::vector<ComplexType> x(hi-lo);
  this->getElementRange(lo, hi, &x[0]);
  for (std::vector<ComplexType>::iterator i = x.begin();
       i != x.end(); ++i) {
    *i = imag(*i);
  }
  this->setElementRange(lo, hi, &x[0]);
  this->ready();
}

// -------------------------------------------------------------
// VectorImplementation::p_conjugate
// -------------------------------------------------------------
void
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

