/**
 * @file   network_partition.cpp
 * @author William A. Perkins
 * @date   2013-08-01 12:30:42 d3g096
 * 
 * @brief  A test of network partitioning
 * 
 * 
 */

#include <iostream>
#include <ga++.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "gridpack/component/base_component.hpp"
#include "base_network.hpp"

// -------------------------------------------------------------
//  class BogusBus
// -------------------------------------------------------------
class BogusBus 
  : public gridpack::component::BaseBusComponent {
public:

  /// Default constructor.
  BogusBus(void)
    : gridpack::component::BaseBusComponent()
  {}

  /// Destructor
  ~BogusBus(void)
  {}

private:

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<BaseBusComponent>(*this);
  }  
};

BOOST_CLASS_EXPORT(BogusBus);

// -------------------------------------------------------------
//  class BogusBranch
// -------------------------------------------------------------
class BogusBranch 
  : public gridpack::component::BaseBranchComponent 
{
public:

  /// Default constructor.
  BogusBranch(void) 
    : gridpack::component::BaseBranchComponent()
  {}

  /// Destructor
  ~BogusBranch(void)
  {}

private:

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<BaseBranchComponent>(*this);
  }  
};

BOOST_CLASS_EXPORT(BogusBranch);

// -------------------------------------------------------------
//  class BogusNetwork
// -------------------------------------------------------------
/// A simple network
/**
 * This builds a simple linear network on the root process.
 * 
 */
class BogusNetwork 
  : public gridpack::network::BaseNetwork<BogusBus, BogusBranch>
{
public:

  /// Default constructor.
  BogusNetwork(const gridpack::parallel::Communicator& comm, const int& local_size)
    : gridpack::network::BaseNetwork<BogusBus, BogusBranch>(comm)
  {
    int global_buses(local_size*this->processor_size());
    int global_branches(global_buses - 1);

    // put all components on process 0
    if (this->processor_rank() == 0) {
      for (int busidx = 0; busidx < global_buses; ++busidx) {
        this->addBus(busidx);
        this->setGlobalBusIndex(busidx, busidx);
      }
      for (int branchidx = 0; branchidx < global_branches; ++branchidx) {
        int bus1(branchidx), bus2(bus1+1);
        this->addBranch(bus1, bus2);
        this->setGlobalBranchIndex(branchidx, branchidx);
        this->setGlobalBusIndex1(branchidx, bus1);
        this->setGlobalBusIndex2(branchidx, bus2);
      }
    }
  }

  /// Destructor
  ~BogusNetwork(void)
  {}

  /// Print the global indexes of the locally owned buses
  void print_bus_ids()
  {
    for (int p = 0; p < this->processor_size(); ++p) {
      if (p == this->processor_rank()) {
        std::cout << p << ": buses: ";
        for (int b = 0; b < this->numBuses(); ++b) {
          int busidx(this->getGlobalBusIndex(b));
          std::cout << busidx;
          if (!this->getActiveBus(b))
            std::cout << '*';
          std::cout << ", ";
        }
        std::cout << std::endl;
        std::cout.flush();
      }
      (this->communicator()).barrier();
    }
  }

  /// Print the global indexes of the locally owned branches
  void print_branch_ids()
  {
    for (int p = 0; p < this->processor_size(); ++p) {
      if (p == this->processor_rank()) {
        std::cout << p << ": branches: ";
        for (int b = 0; b < this->numBranches(); ++b) {
          int branchidx(this->getGlobalBranchIndex(b));
          std::cout << branchidx;
          if (!this->getActiveBranch(b))
            std::cout << '*';
          std::cout << ", ";
        }
        std::cout << std::endl;
        std::cout.flush();
      }
      (this->communicator()).barrier();
    }
  }
};

BOOST_AUTO_TEST_SUITE ( network ) 

BOOST_AUTO_TEST_CASE ( bus_data_serialization )
{
  gridpack::parallel::Communicator world;
  gridpack::network::BusData<BogusBus> bus1;
  
  int an_int_value(14);
  char key[] = "AValue";


  if (world.rank() == 0) {
    bus1.p_originalBusIndex = an_int_value;
    bus1.p_globalBusIndex = an_int_value;
    bus1.p_branchNeighbors.push_back(an_int_value);
    bus1.p_refFlag = true;
    bus1.p_bus->setReferenceBus(true);

    boost::shared_ptr<gridpack::component::DataCollection> 
      busdata(new gridpack::component::DataCollection());
    busdata->addValue(key, an_int_value);
    bus1.p_data = busdata;
  }

  broadcast(world, &bus1, 1, 0);

  BOOST_CHECK_EQUAL(bus1.p_originalBusIndex, an_int_value);
  BOOST_CHECK_EQUAL(bus1.p_globalBusIndex, an_int_value);
  BOOST_CHECK_EQUAL(bus1.p_branchNeighbors.back(), an_int_value);
  BOOST_CHECK(bus1.p_refFlag);
  BOOST_CHECK(bus1.p_data);
  BOOST_CHECK(bus1.p_bus->getReferenceBus());

  int junk;
  bus1.p_data->getValue(key, &junk);
  BOOST_CHECK_EQUAL(junk, an_int_value);
}

BOOST_AUTO_TEST_CASE (bus_data_shuffle)
{
  gridpack::parallel::Communicator world;
  const int local_size(3);
  typedef gridpack::network::BusData<BogusBus> BogusBusData;
  typedef boost::shared_ptr<BogusBusData> BogusBusDataPtr;
  typedef std::vector< BogusBusDataPtr > BogusBusVector;
  char key[] = "AValue";
  
  BogusBusVector buses;
  std::vector<int> dest;
  if (world.rank() == 0) {
    dest.resize(local_size*world.size());
    for (int i = 0; i < local_size*world.size(); ++i) {
      BogusBusDataPtr bus1(new BogusBusData());

      bus1->p_originalBusIndex = i;
      bus1->p_globalBusIndex = i;
      bus1->p_branchNeighbors.push_back(i);
      bus1->p_refFlag = true;
      bus1->p_bus->setReferenceBus(true);
      
      boost::shared_ptr<gridpack::component::DataCollection> 
        busdata(new gridpack::component::DataCollection());
      busdata->addValue(key, i);
      bus1->p_data = busdata;
      buses.push_back(bus1);
      dest[i] = i % world.size();
    }
  }

  Shuffler<BogusBusDataPtr> shuffler;
  shuffler(world, buses, dest);

}  

BOOST_AUTO_TEST_CASE (branch_data_shuffle)
{
  gridpack::parallel::Communicator world;
  const int local_size(3);
  typedef gridpack::network::BranchData<BogusBranch> BogusBranchData;
  typedef boost::shared_ptr<BogusBranchData> BogusBranchDataPtr;
  typedef std::vector< BogusBranchDataPtr > BogusBranchVector;
  char key[] = "AValue";
  
  BogusBranchVector branches;
  std::vector<int> dest;
  if (world.rank() == 0) {
    dest.resize(local_size*world.size());
    for (int i = 0; i < local_size*world.size(); ++i) {
      BogusBranchDataPtr branch1(new BogusBranchData());

      branch1->p_globalBranchIndex = i;
      branch1->p_originalBusIndex1 = i;
      branch1->p_originalBusIndex2 = i;
      branch1->p_globalBusIndex1 = i;
      branch1->p_globalBusIndex2 = i;

      boost::shared_ptr<gridpack::component::DataCollection> 
        branchdata(new gridpack::component::DataCollection());
      branchdata->addValue(key, i);
      branch1->p_data = branchdata;
      branches.push_back(branch1);
      dest[i] = i % world.size();
    }
  }

  Shuffler<BogusBranchDataPtr> shuffler;
  shuffler(world, branches, dest);

}  

BOOST_AUTO_TEST_CASE (branch_data_serialization )
{
  gridpack::parallel::Communicator world;
  gridpack::network::BranchData<BogusBranch> branch1;
  
  int an_int_value(14);
  char key[] = "AValue";

  boost::shared_ptr<gridpack::component::DataCollection> 
    branchdata(new gridpack::component::DataCollection());
  branchdata->addValue(key, an_int_value);

  if (world.rank() == 0) {
    branch1.p_globalBranchIndex = an_int_value;
    branch1.p_globalBusIndex1 = an_int_value;
    branch1.p_globalBusIndex2 = an_int_value;
    branch1.p_originalBusIndex1 = an_int_value;
    branch1.p_originalBusIndex2 = an_int_value;
    branch1.p_activeBranch = false;
    branch1.p_data = branchdata;
  }

  broadcast(world, &branch1, 1, 0);

  BOOST_CHECK_EQUAL(branch1.p_globalBranchIndex, an_int_value);
  BOOST_CHECK_EQUAL(branch1.p_globalBusIndex1, an_int_value);
  BOOST_CHECK_EQUAL(branch1.p_globalBusIndex2, an_int_value);
  BOOST_CHECK_EQUAL(branch1.p_originalBusIndex1, an_int_value);
  BOOST_CHECK_EQUAL(branch1.p_originalBusIndex2, an_int_value);
  BOOST_CHECK(!branch1.p_activeBranch);
  BOOST_CHECK(branch1.p_data);

  int junk;
  branch1.p_data->getValue(key, &junk);
  BOOST_CHECK_EQUAL(junk, an_int_value);

}

BOOST_AUTO_TEST_CASE ( partition )
{
  gridpack::parallel::Communicator world;
  static const int local_size(3);
  BogusNetwork net(world, local_size);

  int allbuses(net.totalBuses());
  int locbuses(net.numBuses());
  BOOST_CHECK_EQUAL(allbuses, local_size*world.size());
  if (world.rank() == 0) {
    BOOST_CHECK_EQUAL(locbuses, allbuses);
  } else {
    BOOST_CHECK_EQUAL(locbuses, 0);
  }

  net.print_bus_ids();
  net.print_branch_ids();
  net.write_graph("network-before.dot");

  net.partition();

  net.print_bus_ids();
  net.print_branch_ids();
  net.write_graph("network-after.dot");
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



