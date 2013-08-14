/**
 * @file   pf_network_test.cpp
 * @author William A. Perkins
 * @date   2013-08-14 12:24:49 d3g096
 * 
 * @brief  Unit tests for powerflow network and component types
 * 
 * 
 */

#include <iostream>
#include <ga++.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/scoped_ptr.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "gridpack/component/base_component.hpp"
#include "pf_components.hpp"

BOOST_AUTO_TEST_SUITE ( pf_network ) 

BOOST_AUTO_TEST_CASE ( bus_serialization )
{
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::powerflow::PFBus> bus1;

  if (world.rank() == 0) {
    bus1.reset(new gridpack::powerflow::PFBus);
    // FIXME: need initialization
  }
  broadcast(world, bus1, 0);

  // FIXME: need to test contents
}

BOOST_AUTO_TEST_CASE ( branch_serialization )
{
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::powerflow::PFBranch> branch1;

  if (world.rank() == 0) {
    branch1.reset(new gridpack::powerflow::PFBranch);
    // FIXME: need initialization
  }
  broadcast(world, branch1, 0);

  // FIXME: need to test contents
}

BOOST_AUTO_TEST_CASE ( partition )
{
  gridpack::parallel::Communicator world;
  static const int local_size(3);
  boost::scoped_ptr<gridpack::powerflow::PFNetwork> 
    net(new gridpack::powerflow::PFNetwork(world));

  int global_buses(local_size*net->processor_size());
  int global_branches(global_buses - 1);

  // Build a simple linear network and put all components on process 0

  if (net->processor_rank() == 0) {
    for (int busidx = 0; busidx < global_buses; ++busidx) {
      net->addBus(busidx);
      net->setGlobalBusIndex(busidx, busidx);
    }
    for (int branchidx = 0; branchidx < global_branches; ++branchidx) {
      int bus1(branchidx), bus2(bus1+1);
      net->addBranch(bus1, bus2);
      net->setGlobalBranchIndex(branchidx, branchidx);
      net->setGlobalBusIndex1(branchidx, bus1);
      net->setGlobalBusIndex2(branchidx, bus2);
    }
  }

  int allbuses(net->totalBuses());
  int locbuses(net->numBuses());
  BOOST_CHECK_EQUAL(allbuses, global_buses);
  if (world.rank() == 0) {
    BOOST_CHECK_EQUAL(locbuses, allbuses);
  } else {
    BOOST_CHECK_EQUAL(locbuses, 0);
  }

  net->write_graph("network-before.dot");

  net->partition();

  allbuses = net->totalBuses();
  locbuses = net->numBuses();
  BOOST_CHECK_EQUAL(allbuses, global_buses);
  BOOST_CHECK(locbuses > 0);

  net->write_graph("network-after.dot");
}

BOOST_AUTO_TEST_SUITE_END( )

// -------------------------------------------------------------
// init_function
// -------------------------------------------------------------
bool init_function()
{
  return true;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  GA_Initialize();
  MA_init(MT_C_CHAR, 1024*1024, 1024*1024);
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  GA::Terminate();
}


