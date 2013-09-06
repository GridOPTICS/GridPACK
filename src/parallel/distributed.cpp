// -------------------------------------------------------------
/**
 * @file   distributed.cpp
 * @author William A. Perkins
 * @date   2013-09-06 14:57:42 d3g096
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
//  class DistributedInterface
// -------------------------------------------------------------

// -------------------------------------------------------------
// DistributedInterface:: constructors / destructor
// -------------------------------------------------------------
DistributedInterface::DistributedInterface(void)
{
  // empty
}

DistributedInterface::DistributedInterface(const DistributedInterface& old)
{
  // empty
}

DistributedInterface::~DistributedInterface(void)
{
}

// -------------------------------------------------------------
//  class Distributed
// -------------------------------------------------------------

// -------------------------------------------------------------
// Distributed:: constructors / destructor
// -------------------------------------------------------------
Distributed::Distributed(const Communicator& comm)
  : DistributedInterface(), communicator_(comm)
{
  // empty
}

Distributed::Distributed(const Distributed& old)
  : DistributedInterface(old), communicator_(old.communicator_)
{
  // empty
}

Distributed::~Distributed(void)
{
  // empty
}

// -------------------------------------------------------------
// Distributed::communicator
// -------------------------------------------------------------
const Communicator& 
Distributed::communicator(void) const
{
  return communicator_;
}

// -------------------------------------------------------------
//  class WrappedDistributed
// -------------------------------------------------------------

// -------------------------------------------------------------
// WrappedDistributed:: constructors / destructor
// -------------------------------------------------------------
WrappedDistributed::WrappedDistributed(Distributed *d)
  : DistributedInterface(), p_distributed(d)
{
  
}

WrappedDistributed::WrappedDistributed(const WrappedDistributed& old)
  : p_distributed(old.p_distributed)
{
  
}

WrappedDistributed::~WrappedDistributed(void)
{
}

// -------------------------------------------------------------
// WrappedDistributed::communicator
// -------------------------------------------------------------
const Communicator& 
WrappedDistributed::communicator(void) const
{
  return p_distributed->communicator();
}





} // namespace parallel
} // namespace gridpack

