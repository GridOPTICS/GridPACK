// -------------------------------------------------------------
/**
 * @file   test_network.cpp
 * @author Bruce Palmer, William Perkins
 * @date   2013-07-18 12:09:54 d3g096
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
  gridpack::parallel::Communicator world;
  gridpack::network::BaseNetwork<int,int> network(world);
  gridpack::network::BusData<int> bus1;
  gridpack::network::BusData<int> bus2;
  bus1 = bus2;
  gridpack::network::BranchData<int> branch1;
  gridpack::network::BranchData<int> branch2;
  branch1 = branch2;
}

/**
 * Default destructor.
 */
gridpack::TestNetwork::~TestNetwork(void)
{
}
