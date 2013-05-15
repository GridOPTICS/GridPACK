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

#include "gridpack/parallel/distribution.hpp"
#include "gridpack/base_component/base_network_component.hpp"

/**
 * Simple constructor
 */
BaseNetworkComponent::BaseNetworkComponent()
{
  p_network = static_cast<*BaseNetwork>0;
  p_idx = -1;
}

/**
 * Constructor
 * @param network: pointer to network the component is associated with
 * @param idx: local bus or branch index that network is associated with
 */
BaseNetworkComponent::BaseNetworkComponent(BaseNetwork *network, int idx)
{
  p_network = network;
  p_idx = idx;
}

/**
 * Constructor without local network index
 * @param network: pointer to network the component is associated with
 */
BaseNetworkComponent::BaseNetworkComponent(BaseNetwork *network)
{
  p_network = network;
  p_idx = -1;
}

/**
 * Destructor
 */
BaseNetworkComponent::~BaseNetworkComponent(void)
{
  p_network = static_cast<stlplus::smart_ptr<BaseNetwork> >0;
  p_idx = -1;
}

/**
 * Set the network associated with the component
 */
void BaseNetworkComponent::setNetwork(stlplus::smart_ptr<BaseNetwork> network)
{
  p_network = network;
}

/**
 * Set the value of the local network index the component is associated with
 * @param idx: value of local network index
 */
void BaseNetworkComponent::setIndex(int idx)
{
  p_idx = idx;
}

/**
 * Return the size of the component for use in packing and
 * unpacking routines. This might not be needed, but throw
 * it in for now.
 * @return: size of network component
 */
int BaseNetworkComponent::size(void)
{
  return sizeof(*this);
}

