// -------------------------------------------------------------
/**
 * @file   distributed.cpp
 * @author William A. Perkins
 * @date   2013-05-07 07:48:09 d3g096
 * 
 * @brief  Implementation of gridpack::parallel::Distributed
 * 
 * 
 */
// -------------------------------------------------------------

#include "distributed.hpp"

namespace gridpack {
namespace parallel {

// -------------------------------------------------------------
//  class Distributed
// -------------------------------------------------------------

// -------------------------------------------------------------
// Distributed:: constructors / destructor
// -------------------------------------------------------------
Distributed::Distributed(const Communicator& comm)
  : communicator_(comm)
{
  // empty
}

Distributed::Distributed(const Distributed& old)
  : communicator_(old.communicator_)
{
  // empty
}

Distributed::~Distributed(void)
{
  // empty
}

// -------------------------------------------------------------
// Distributed
// -------------------------------------------------------------


} // namespace parallel
} // namespace gridpack

