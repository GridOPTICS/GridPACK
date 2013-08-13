// -------------------------------------------------------------
/**
 * @file   vector_implementation.cpp
 * @author William A. Perkins
 * @date   2013-08-13 12:20:06 d3g096
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
// VectorImplementation::p_set_element_range
// -------------------------------------------------------------
void
VectorImplementation::p_set_element_range(const int& lo, const int& hi, ComplexType *x)
{
  std::vector<int> i;
  i.reserve(hi-lo);
  std::copy(boost::counting_iterator<int>(lo),
            boost::counting_iterator<int>(hi),
            std::back_inserter(i));
  this->p_set_elements(i.size(), &i[0], x);
}

// -------------------------------------------------------------
// VectorImplementation::p_get_element_range
// -------------------------------------------------------------
void
VectorImplementation::p_get_element_range(const int& lo, const int& hi, ComplexType *x) const
{
  std::vector<int> i;
  i.reserve(hi-lo);
  std::copy(boost::counting_iterator<int>(lo),
            boost::counting_iterator<int>(hi),
            std::back_inserter(i));
  this->p_get_elements(i.size(), &i[0], x);
}



} // namespace math
} // namespace gridpack

