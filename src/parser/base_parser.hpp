/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * GOSS_parser.hpp
 *
 *  Created on: December 30, 2014
 *      Author: Kevin Glass, Bruce Palmer
 */

#ifndef BASEPARSER_HPP_
#define BASEPARSER_HPP_

#define OLD_MAP

#include <iostream>
#include <string>
#include <vector>
#include <map>
#ifndef OLD_MAP
#include <boost/unordered_map.hpp>
#endif


#include "gridpack/component/base_component.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/exception.hpp"
#include "gridpack/network/base_network.hpp"

namespace gridpack {
namespace parser {

template <class _network>
class BaseParser
{
  public:

    /// Constructor
    explicit BaseParser() : p_configExists(false) {}

    /**
     * Destructor
     */
    virtual ~BaseParser(){}

    void createNetwork(std::vector<boost::shared_ptr<component::DataCollection> >
        &busData, std::vector<boost::shared_ptr<component::DataCollection> >
        &branchData)
    {
      p_timer = gridpack::utility::CoarseTimer::instance();
      int t_create = p_timer->createCategory("Parser:createNetwork");
      p_timer->start(t_create);
      int me(p_network->communicator().rank());
      int nprocs(p_network->communicator().size());
      int i;
      // Exchange information on number of buses and branches on each
      // processor
      std::vector<int> sbus(nprocs);
      std::vector<int> sbranch(nprocs);
      std::vector<int> nbus(nprocs);
      std::vector<int> nbranch(nprocs);
      for (i=0; i<nprocs; i++) {
        sbus[i] = 0;
        sbranch[i] = 0;
      }
      sbus[me] = busData.size();
      sbranch[me] = branchData.size();
      MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
      int ierr;
      ierr = MPI_Allreduce(&sbus[0],&nbus[0],nprocs,MPI_INT,MPI_SUM,comm);
      ierr = MPI_Allreduce(&sbranch[0],&nbranch[0],nprocs,MPI_INT,MPI_SUM,comm);
      // evaluate offsets for buses and branches
      std::vector<int> offset_bus(nprocs);
      std::vector<int> offset_branch(nprocs);
      offset_bus[0] = 0;
      offset_branch[0] = 0;
      for (i=1; i<nprocs; i++) {
        offset_bus[i] = offset_bus[i-1]+nbus[i-1];
        offset_branch[i] = offset_branch[i-1]+nbranch[i-1];
      }

      int numBus = busData.size();
      for (i=0; i<numBus; i++) {
        int idx;
        busData[i]->getValue(BUS_NUMBER,&idx);
        p_network->addBus(idx);
        p_network->setGlobalBusIndex(i,i+offset_bus[me]);
        *(p_network->getBusData(i)) = *(busData[i]);
        p_network->getBusData(i)->addValue(CASE_ID,p_case_id);
        p_network->getBusData(i)->addValue(CASE_SBASE,p_case_sbase);
      }
      int numBranch = branchData.size();
      for (i=0; i<numBranch; i++) {
        int idx1, idx2;
        branchData[i]->getValue(BRANCH_FROMBUS,&idx1);
        branchData[i]->getValue(BRANCH_TOBUS,&idx2);
        p_network->addBranch(idx1, idx2);
        p_network->setGlobalBranchIndex(i,i+offset_branch[me]);
#ifdef OLD_MAP
        std::map<int, int>::iterator it;
#else
        boost::unordered_map<int, int>::iterator it;
#endif
        *(p_network->getBranchData(i)) = *(branchData[i]);
        p_network->getBranchData(i)->addValue(CASE_ID,p_case_id);
        p_network->getBranchData(i)->addValue(CASE_SBASE,p_case_sbase);
      }
      p_configExists = true;
#if 0
      // debug
      printf("Number of buses: %d\n",numBus);
      for (i=0; i<numBus; i++) {
        printf("Dumping bus: %d\n",i);
        p_network->getBusData(i)->dump();
      }
      printf("Number of branches: %d\n",numBranch);
      for (i=0; i<numBranch; i++) {
        printf("Dumping branch: %d\n",i);
        p_network->getBranchData(i)->dump();
      }
#endif
      busData.clear();
      branchData.clear();
      p_timer->stop(t_create);
    }

  protected:

    /* ************************************************************************
     **************************************************************************
     ***** PROTECTED SCOPE
     **************************************************************************
     *********************************************************************** */
    /**
     * Assign network to internal network pointer variable
     */
    void setNetwork(boost::shared_ptr<_network> network)
    {
      p_network = network;
      p_network_data = network->getNetworkData();
    }

    /**
     * Set some global parameters
     */
    void setCaseID(int id)
    {
      p_case_id = id;
    }

    void setCaseSBase(double sb)
    {
      p_case_sbase = sb;
    }

    bool configExists(void)
    {
      return p_configExists;
    }

    /* ************************************************************************
     **************************************************************************
     ***** OBJECT DATA
     **************************************************************************
     *********************************************************************** */
    boost::shared_ptr<_network> p_network;

    int                      p_case_id;
    double                   p_case_sbase;
    gridpack::utility::CoarseTimer *p_timer;
    bool                     p_configExists;
    /**
     * Data collection object associated with network as a whole
     */
    boost::shared_ptr<gridpack::component::DataCollection> p_network_data;
}; /* end of GOSS_parser */

} /* namespace parser */
} /* namespace gridpack */

#endif /* BASEPARSER_HPP_ */
