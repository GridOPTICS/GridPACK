// -------------------------------------------------------------
/**
 * @file   base_network_component.cpp
 * @author Bruce Palmer
 * @date   May 13, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/parallel/distributed.hpp"
#include "gridpack/base_component/base_network_component.hpp"

/**
 * Simple constructor
 */
BaseNetworkComponent::BaseNetworkComponent()
{
#if 0
  p_network = static_cast<*BaseNetwork>0;
#endif
  p_idx = -1;
}

/**
 * Constructor
 * @param network: pointer to network the component is associated with
 * @param idx: local bus or branch index that network is associated with
 */
#if 0
BaseNetworkComponent::BaseNetworkComponent(BaseNetwork *network, int idx)
{
  p_network = network;
  p_idx = idx;
}
#endif

/**
 * Constructor without local network index
 * @param network: pointer to network the component is associated with
 */
#if 0
BaseNetworkComponent::BaseNetworkComponent(BaseNetwork *network)
{
  p_network = network;
  p_idx = -1;
}
#endif

/**
 * Destructor
 */
BaseNetworkComponent::~BaseNetworkComponent(void)
{
#if 0
  p_network = static_cast<boost::smart_ptr<BaseNetwork> >0;
#endif
  p_idx = -1;
}

/**
 * Set the network associated with the component
 */
#if 0
void BaseNetworkComponent::setNetwork(boost::smart_ptr<BaseNetwork> network)
{
  p_network = network;
}
#endif

/**
 * Set the value of the local network index the component is associated with
 * @param idx: value of local network index
 */
void BaseNetworkComponent::setIndex(int idx)
{
  p_idx = idx;
}

/**
 * Get the value of the local network index the component is associated with
 * @return: value of local network index
 */
int getIndex()
{
  return p_idx;
}

/**
 * Return the size of the component for use in packing and
 * unpacking routines. This might not be needed, but throw
 * it in for now.
 * @return: size of network component
 */
int BaseNetworkComponent::size(void) const
{
  return sizeof(*this);
}

