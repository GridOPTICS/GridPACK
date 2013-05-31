// -------------------------------------------------------------
/**
 * @file   test_network.cpp
 * @author Bruce Palmer, William Perkins
 * @date   April 3, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/network/test_network.hpp"
#include "gridpack/network/base_network.hpp"

// -------------------------------------------------------------
//  class TestNetwork:
//  Trivial class that instantiates network so we can check for compile time
//  errors in BaseNetwork class.
// -------------------------------------------------------------

/**
 * Constructor.
 */
gridpack::TestNetwork::TestNetwork(void)
{
  gridpack::network::BaseNetwork<int, int> network;
};

/**
 * Default destructor.
 */
gridpack::TestNetwork::~TestNetwork(void)
{
};
