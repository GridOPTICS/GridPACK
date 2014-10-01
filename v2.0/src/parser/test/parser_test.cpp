/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   parser_test.cpp
 * @author Bruce Palmer
 * @date   July 29, 2014
 * 
 * @brief  A simple test of PTI23_parser test and HashDistribution
 *         functionality
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------

#include <iostream>
#include <ga.h>
#include "gridpack/parallel/communicator.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/parser/hash_distr.hpp"

class TestBus
  : public gridpack::component::BaseBusComponent {
  public:

  TestBus(void) {
  }

  ~TestBus(void) {
  }
};

class TestBranch
  : public gridpack::component::BaseBranchComponent {
  public:

  TestBranch(void) {
  }

  ~TestBranch(void) {
  }
};

typedef gridpack::network::BaseNetwork<TestBus, TestBranch> TestNetwork;

struct bus_data {int idx;};
struct branch_data {int idx1; int idx2;};

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  GA_Initialize();
  // Create an artificial scope so that all objects call their destructors
  // before GA_Terminate is called
  if (1) {
    gridpack::parallel::Communicator world;
    MPI_Comm comm = static_cast<MPI_Comm>(world);
    boost::shared_ptr<TestNetwork> network(new TestNetwork(world));
    // Create parser and parse data
    gridpack::parser::PTI23_parser<TestNetwork> parser(network);
    parser.parse("parser_data.raw");
    // Distribute data across processors
    network->partition();

    // Check to see if contents of DataCollection objects is correct
    int nbus = network->numBuses();
    int i, idx;
    int schk, rchk;
    double rval;
    schk = 0;
    for (i=0; i<nbus; i++) {
      network->getBusData(i)->getValue(LOAD_PL,&rval);
      idx = network->getGlobalBusIndex(i);
      if (rval != static_cast<double>(idx%100)) schk = 1;
      network->getBusData(i)->getValue(LOAD_QL,&rval);
      if (rval != static_cast<double>((idx+1)%100)) schk = 1;
    }

    int nbranch = network->numBranches();
    for (i=0; i<nbranch; i++) {
      network->getBranchData(i)->getValue(BRANCH_RATING_A,&rval,0);
      idx = network->getGlobalBranchIndex(i);
      if (rval != static_cast<double>(idx%100)) schk = 1;
    }
    MPI_Allreduce(&schk,&rchk,1,MPI_INT,MPI_SUM,comm);
    if (rchk == 0 && world.rank() == 0) {
      printf("\nParsing of test configuration is ok\n");
    } else if (world.rank() == 0) {
      printf("\nError in parsing of test configuration\n");
    }

    // Check to see if hash distribution functionality works
    gridpack::hash_distr::HashDistribution<TestNetwork,bus_data,branch_data>
      hashMap(network);
    // create data on process 0 to go to all buses in the network
    int bus_total = network->totalBuses();
    std::vector<bus_data> busData;
    std::vector<int> keys;
    if (world.rank() == 0) {
      for (i=0; i<bus_total; i++) {
        bus_data bus;
        bus.idx = i+1;
        busData.push_back(bus);
        keys.push_back(i+1);
      }
    }
    hashMap.distributeBusValues(keys,busData);
    schk = 0;
    for (i=0; i<keys.size(); i++) {
      idx = keys[i];
      if (network->getGlobalBusIndex(idx)+1 != busData[i].idx) schk = 1;
    }
#if 1
    // create data on process 0 to go to all branches in the network
    int branch_total = network->totalBranches();
    int dim = static_cast<int>(sqrt(static_cast<double>(bus_total))+0.0000001);
    //int dim = sqrt(bus_total);
    std::vector<branch_data> branchData;
    keys.clear();
    std::vector<std::pair<int,int> > pairs;
    int ix, iy, idx1, idx2;
    if (world.rank() == 0) {
      for (i=0; i<branch_total; i++) {
        idx = i;
        branch_data branch;
        branch.idx1 = idx+1;
        branch.idx2 = idx+2;
        branchData.push_back(branch);
        if (idx < (dim-1)*dim) {
          ix = idx%(dim-1);
          iy = (idx-ix)/(dim-1);
          ix++;
          idx1 = iy*dim+ix;
          idx2 = iy*dim+ix+1;
        } else {
          idx -= (dim-1)*dim;
          ix = idx%dim;
          iy = (idx-ix)/dim;
          iy++;
          idx1 = (iy-1)*dim+ix+1;
          idx2 = iy*dim+ix+1;
        }
        pairs.push_back(std::pair<int,int>(idx1,idx2));
      }
    }
    hashMap.distributeBranchValues(pairs,keys,branchData);
    for (i=0; i<keys.size(); i++) {
      idx = keys[i];
      idx = network->getGlobalBranchIndex(idx);
      if (idx+1 != branchData[i].idx1 || idx+2 != branchData[i].idx2) schk = 1;
    }
#endif
    MPI_Allreduce(&schk,&rchk,1,MPI_INT,MPI_SUM,comm);
    if (rchk == 0 && world.rank() == 0) {
      printf("\nHash distribution is ok\n");
    } else if (world.rank() == 0) {
      printf("\nError in hash distribution\n");
    }
  }

  GA_Terminate();
  return 0;
}

