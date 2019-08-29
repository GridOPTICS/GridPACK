/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_factory.hpp
 * @author Bruce Palmer
 * @date   2013-08-08 10:21:10 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _base_factory_h_
#define _base_factory_h_

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"

// Base factory class that contains functions that are generic to all
// applications.

namespace gridpack{
namespace factory{

template <class _network>
class BaseFactory {
  public:
    typedef _network NetworkType;
    typedef boost::shared_ptr<NetworkType> NetworkPtr;

    /**
     * Constructor
     * @param network network that factory operates on
     */
    BaseFactory(NetworkPtr network)
      : p_network(network)
    { 
      p_profile = false;
      p_numBuses = p_network->numBuses();
      p_numBranches = p_network->numBranches();
      p_buses = new gridpack::component::BaseBusComponent*[p_numBuses];
      p_branches = new gridpack::component::BaseBranchComponent*[p_numBranches];
      int i;
      for (i=0; i<p_numBuses; i++) {
        p_buses[i] = p_network->getBus(i).get();
      }
      for (i=0; i<p_numBranches; i++) {
        p_branches[i] = p_network->getBranch(i).get();
      }
    }

    /**
     * Destructor
     */
    virtual ~BaseFactory(void)
    {
      delete [] p_buses;
      delete [] p_branches;
    }

    /**
     * Set pointers in each bus and branch component so that it points to
     * connected buses and branches. This routine operates on the generic
     * BaseBusComponent and BaseBranchComponent interfaces. It also sets some
     * indices in MatVecInterface for each component.
     */
    virtual void setComponents(void)
    {
      gridpack::utility::CoarseTimer *timer =
        gridpack::utility::CoarseTimer::instance();
      timer->configTimer(p_profile);
      int t_setc = timer->createCategory("Factory:setComponents");
      timer->start(t_setc);
      int i, j;
      int idx1, idx2;

      // Set pointers for buses at either end of each branch
      int numActiveBranch = 0;
      for (i=0; i<p_numBranches; i++) {
        int branch_idx, bus1_idx, bus2_idx;
        p_network->getBranchEndpoints(i, &idx1, &idx2);
        p_network->getBranch(i)->setBus1(p_network->getBus(idx1));
        p_network->getBranch(i)->setBus2(p_network->getBus(idx2));
        branch_idx = p_network->getGlobalBranchIndex(i); 
        bus1_idx = p_network->getGlobalBusIndex(idx1); 
        bus2_idx = p_network->getGlobalBusIndex(idx2); 
        p_network->getBranch(i)->setMatVecIndices(bus1_idx, bus2_idx); 
        if (p_network->getActiveBranch(i)) numActiveBranch++;
      }

      // Set pointers for branches and buses connected to each bus
      int numActiveBus = 0;
      for (i=0; i<p_numBuses; i++) {
        p_network->getBus(i)->clearBuses();
        std::vector<int> nghbrBus = p_network->getConnectedBuses(i);
        for (j=0; j<nghbrBus.size(); j++) {
          p_network->getBus(i)->addBus(p_network->getBus(nghbrBus[j]));
        }
        p_network->getBus(i)->clearBranches();
        std::vector<int> nghbrBranch = p_network->getConnectedBranches(i);
        for (j=0; j<nghbrBranch.size(); j++) {
          p_network->getBus(i)->addBranch(p_network->getBranch(nghbrBranch[j]));
        }
        int bus_idx;
        bus_idx = p_network->getGlobalBusIndex(i); 
        p_network->getBus(i)->setMatVecIndex(bus_idx);
        if (p_network->getActiveBus(i)) numActiveBus++;
      }
      
      // Set reference bus
      int idx = p_network->getReferenceBus();
      if (idx != -1) {
        p_network->getBus(idx)->setReferenceBus(true);
      }

      // Set bus and branch indices
      for (i=0; i<p_numBuses; i++) {
        p_network->getBus(i)->setOriginalIndex(p_network->getOriginalBusIndex(i));
        p_network->getBus(i)->setGlobalIndex(p_network->getGlobalBusIndex(i));
      }
      for (i=0; i<p_numBranches; i++) {
        gridpack::component::BaseBusComponent *bus1 =
          dynamic_cast<gridpack::component::BaseBusComponent*>
          (p_network->getBranch(i)->getBus1().get());
        gridpack::component::BaseBusComponent *bus2 =
          dynamic_cast<gridpack::component::BaseBusComponent*>
          (p_network->getBranch(i)->getBus2().get());
        p_network->getBranch(i)->setGlobalIndex(
            p_network->getGlobalBranchIndex(i));
        p_network->getBranch(i)->setBus1OriginalIndex(
            bus1->getOriginalIndex());
        p_network->getBranch(i)->setBus2OriginalIndex(
            bus2->getOriginalIndex());
        p_network->getBranch(i)->setBus1GlobalIndex(
            bus1->getGlobalIndex());
        p_network->getBranch(i)->setBus2GlobalIndex(
            bus2->getGlobalIndex());
      }

      // Come up with a set of global indices for each component so that the buses
      // and branches are consecutively numbered on each component and the indices
      // for all active components are unique. These are used in the mapper routine.
      // First get the number of active buses and branches on each process and
      // broadcast this to all other processes.
      int grp = p_network->communicator().getGroup();
      int nprocs = GA_Pgroup_nnodes(grp);
      int me = GA_Pgroup_nodeid(grp);
      int *activeBus = new int[nprocs];
      int *activeBranch = new int[nprocs];

      for (i=0; i<nprocs; i++) {
        activeBus[i] = 0;
        activeBranch[i] = 0;
      }

      activeBus[me] = numActiveBus;
      activeBranch[me] = numActiveBranch;
      char cplus[2];
      strcpy(cplus,"+");
      GA_Pgroup_igop(grp,activeBus,nprocs,cplus);
      GA_Pgroup_igop(grp,activeBranch,nprocs,cplus);

      // Create indices for buses. Start by creating a global array with an entry
      // for each bus
      int one = 1;
      int ntot = 0;
      int offset = 0;
      for (i=0; i<nprocs; i++) {
        ntot += activeBus[i];
        if (i<me) offset += activeBus[i];
      }
      int g_bus = GA_Create_handle();
      GA_Set_data(g_bus, one, &ntot, C_INT);
      GA_Set_pgroup(g_bus, grp);
      if (!GA_Allocate(g_bus)) {
        char buf[256];
        sprintf(buf,"BaseFactory::setComponents: Unable to allocate distributed array for indices");
        printf("%s",buf);
        throw gridpack::Exception(buf);
      }
      GA_Zero(g_bus);

      // Scatter new index values into the global index locations for the buses
      int *ibus_val = new int[numActiveBus];
      int **ibus_idx = new int*[numActiveBus];
      int icnt = 0;
      for (i=0; i<p_numBuses; i++) {
        if (p_network->getActiveBus(i)) {
          ibus_idx[icnt] = new int;
          *(ibus_idx[icnt]) = p_network->getGlobalBusIndex(i);
          ibus_val[icnt] = offset+icnt;
          icnt++;
        }
      }
      NGA_Scatter(g_bus,ibus_val,ibus_idx,numActiveBus);
      GA_Pgroup_sync(grp);
      for (i=0; i<numActiveBus; i++) {
        delete ibus_idx[i];
      }
      delete [] ibus_idx;
      delete [] ibus_val;

      // Now gather values for both active and inactive buses
      ibus_val = new int[p_numBuses];
      ibus_idx = new int*[p_numBuses];
      for (i=0; i<p_numBuses; i++) {
        ibus_idx[i] = new int;
        *(ibus_idx[i]) = p_network->getGlobalBusIndex(i);
      }
      NGA_Gather(g_bus,ibus_val,ibus_idx,p_numBuses);
      // Assign the MatVecIndex for the bus and clean up arrays
      for (i=0; i<p_numBuses; i++) {
        p_network->getBus(i)->setMatVecIndex(ibus_val[i]);
        delete ibus_idx[i];
      }
      delete [] ibus_idx;
      delete [] ibus_val;
      delete [] activeBus;
      delete [] activeBranch;
      GA_Destroy(g_bus);

      // Finish by assigning MatVecIndices for the branches
      for (i=0; i<p_numBranches; i++) {
        p_network->getBranch(i)->getBus1()->getMatVecIndex(&idx1);
        p_network->getBranch(i)->getBus2()->getMatVecIndex(&idx2);
        p_network->getBranch(i)->setMatVecIndices(idx1,idx2);
      }
      
      // Set internal maps
      p_network->setMap();
      timer->stop(t_setc);
      timer->configTimer(true);
    }


    /**
     * Generic method that invokes the "load" method on all branches and buses
     * to move data from the DataCollection objects on the network into the
     * corresponding buses and branches
     */
    virtual void load(void)
    {
      gridpack::utility::CoarseTimer *timer =
        gridpack::utility::CoarseTimer::instance();
      timer->configTimer(p_profile);
      int t_load = timer->createCategory("Factory:load");
      timer->start(t_load);
      int t_nbus = timer->createCategory("Factory:load:nbus");
      timer->start(t_nbus);
      timer->stop(t_nbus);
      int i;
      int rank = p_network->communicator().rank();

      // Invoke load method on all bus objects
      int t_load1 = timer->createCategory("Factory:load:bus");
      timer->start(t_load1);
      for (i=0; i<p_numBuses; i++) {
        p_network->getBus(i)->setRank(rank);
        p_network->getBus(i)->load(p_network->getBusData(i));
        if (p_network->getBus(i)->getReferenceBus())
          p_network->setReferenceBus(i);
      }
      timer->stop(t_load1);

      // Invoke load method on all branch objects
      int t_load2 = timer->createCategory("Factory:load:branch");
      timer->start(t_load2);
      for (i=0; i<p_numBranches; i++) {
        p_network->getBranch(i)->setRank(rank);
        p_network->getBranch(i)->load(p_network->getBranchData(i));
      }
      timer->stop(t_load2);
      timer->stop(t_load);
      timer->configTimer(true);
    }

    /**
     * Set up the exchange buffers so that they work correctly. This should only
     * be called after the network topology has been specified
     * @param flag if true then have network allocate buffers for exchange,
     * otherwise use external buffers
     */
    virtual void setExchange(bool flag=true)
    {
      gridpack::utility::CoarseTimer *timer =
        gridpack::utility::CoarseTimer::instance();
      timer->configTimer(p_profile);
      int t_setx = timer->createCategory("Factory:setExchange");
      timer->start(t_setx);
      int busXCSize, branchXCSize;
      int nbus, nbranch;

      nbus = p_numBuses;
      nbranch = p_numBranches;

      // Get size of bus and branch exchange buffers from first local bus and branch
      // components. These must be the same for all bus and branch components

      if (nbus > 0) {
        busXCSize = p_network->getBus(0)->getXCBufSize();
      } else {
        busXCSize = 0;
      }
      if (flag) {
        p_network->allocXCBus(busXCSize);
      } else {
        p_network->allocXCBusPointers(busXCSize);
      }

      if (nbranch > 0){
        branchXCSize = p_network->getBranch(0)->getXCBufSize();
      } else {
        branchXCSize = 0;
      }
      if (flag) {
        p_network->allocXCBranch(branchXCSize);
      } else {
        p_network->allocXCBranchPointers(branchXCSize);
      }

      int i;
      if (flag) {
        // Buffers have been allocated in network. Now associate buffers from network
        // back to individual components
        if (busXCSize > 0) {
          for (i=0; i<nbus; i++) {
            p_network->getBus(i)->setXCBuf(p_network->getXCBusBuffer(i));
          }
        }
        if (branchXCSize > 0) {
          for (i=0; i<nbranch; i++) {
            p_network->getBranch(i)->setXCBuf(p_network->getXCBranchBuffer(i));
          }
        }
      } else {
        // Buffers have been allocated in components. Now assign buffers to
        // pointers in network
        void *ptr;
        if (busXCSize > 0) {
          for (i=0; i<nbus; i++) {
            p_network->getBus(i)->getXCBuf(&ptr);
            p_network->setXCBusBuffer(i,ptr);
          }
        }
        if (branchXCSize > 0) {
          for (i=0; i<nbranch; i++) {
            p_network->getBranch(i)->getXCBuf(&ptr);
            p_network->setXCBranchBuffer(i,ptr);
          }
        }
      }
      timer->stop(t_setx);
      timer->configTimer(true);
    }

    /**
     * Set the mode for all BaseComponent objects in the network.
     * @param mode integer representing desired mode
     */
    virtual void setMode(int mode)
    {
      int i;
      for (i=0; i<p_numBuses; i++) {
        p_buses[i]->setMode(mode);
      }
      for (i=0; i<p_numBranches; i++) {
        p_branches[i]->setMode(mode);
      }
    }

    /**
     * Set the mode for all BaseBusComponent objects in the network.
     * @param mode integer representing desired mode
     */
    virtual void setBusMode(int mode)
    {
      int i;
      for (i=0; i<p_numBuses; i++) {
        p_buses[i]->setMode(mode);
      }
    }

    /**
     * Set the mode for all BaseBranchComponent objects in the network.
     * @param mode integer representing desired mode
     */
    virtual void setBranchMode(int mode)
    {
      int i;
      for (i=0; i<p_numBranches; i++) {
        p_branches[i]->setMode(mode);
      }
    }

    /**
     * A convenience function that checks to see if something is true on all
     * processors
     * @param flag boolean flag on each processor
     * @return true if flag is true on all processors, false otherwise
     */
    bool checkTrue(bool flag) {
      int iok;
      if (flag) {
        iok = 1;
      } else {
        iok = 0;
      }
      int grp = p_network->communicator().getGroup();
      int nprocs = GA_Pgroup_nnodes(grp);
      char cplus[2];
      strcpy(cplus,"+");
      GA_Pgroup_igop(grp,&iok,1,cplus);
      if (iok == nprocs) {
        return true;
      } else {
        return false;
      }
    }

    /**
     * A convenience function that checks to see if something is true on at
     * least one processor
     * @param flag boolean flag on each processor
     * @return true if flag is true on at least one processor, false otherwise
     */
    bool checkTrueSomewhere(bool flag) {
      int iok;
      if (flag) {
        iok = 1;
      } else {
        iok = 0;
      }
      int grp = p_network->communicator().getGroup();
      int nprocs = GA_Pgroup_nnodes(grp);
      char cplus[2];
      strcpy(cplus,"+");
      GA_Pgroup_igop(grp,&iok,1,cplus);
      if (iok > 0) {
        return true;
      } else {
        return false;
      }
    }

    /**
     * Save internal state variables of the buses and branches to the
     * associated data collection object for possible use in output or to
     * transfer them to another network
     */
    void saveData(void)
    {
      int i;
      // Save data on buses
      for (i=0; i<p_numBuses; i++) {
        p_buses[i]->saveData(p_network->getBusData(i));
      }
      // Save data on branches
      for (i=0; i<p_numBranches; i++) {
        p_branches[i]->saveData(p_network->getBranchData(i));
      }
    }

    /**
     * return a pointer to the bus list
     * @param buses pointer to a list of base bus component pointers
     * @param nbus number of buses on processor
     */
    void getBusPointers(gridpack::component::BaseBusComponent ***buses,
        int *nbus)
    {
      *buses = p_buses;
      *nbus = p_numBuses;
    }

    /**
     * return a pointer to the branch list
     * @param branches pointer to a list of base branch component pointers
     * @param nbranch number of branches on processor
     */
    void getBranchPointers(gridpack::component::BaseBranchComponent ***branches,
        int *nbranch)
    {
      *branches = p_branches;
      *nbranch = p_numBranches;
    }

    /**
     * Debugging call that will dump contents of DataCollection objects on each
     * bus and branch. This call will attempt to guarantee that output is
     * generated sequentially from each processor but this will probably fail
     * because of IO buffering
     */
    void dumpData(void)
    {
      int i, iproc;
      // Write out data to buses
      int nprocs = p_network->communicator().size();
      int me = p_network->communicator().rank();
      for (iproc = 0; iproc<nprocs; iproc++) {
        if (iproc == me) {
          for (i=0; i<p_numBuses; i++) {
            printf("p[%d] Printing data for bus %d\n",me,i);
            p_network->getBusData(i)->dump();
          }
        }
        GA_Sync();
      }
      for (iproc = 0; iproc<nprocs; iproc++) {
        if (iproc == me) {
          for (i=0; i<p_numBranches; i++) {
            printf("p[%d] Printing data for branch %d\n",me,i);
            p_network->getBranchData(i)->dump();
          }
        }
        GA_Sync();
      }
    }

  protected:

    NetworkPtr p_network;

    bool p_profile;

    int p_numBuses;

    int p_numBranches;

    gridpack::component::BaseBusComponent **p_buses;

    gridpack::component::BaseBranchComponent **p_branches;
};

}    // factory
}    // gridpack
#endif
