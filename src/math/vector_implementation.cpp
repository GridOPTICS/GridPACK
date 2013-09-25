// -------------------------------------------------------------
/**
 * @file   vector_implementation.cpp
 * @author William A. Perkins
 * @date   2013-09-25 07:00:13 d3g096
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



} // namespace math
} // namespace gridpack

