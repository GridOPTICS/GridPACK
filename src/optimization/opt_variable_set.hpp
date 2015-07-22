/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   opt_variable_set.hpp
 * @author 
 * @date   
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _opt_variable_set_hpp_
#define _opt_variable_set_hpp_

#include <vector>
#include <macdecls.h>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/optimization_ifc.hpp"
#include <ga.h>


// Base optimization class that contains functions that are generic to all
// applications.

namespace gridpack{
namespace optimization{

template <class _network>
class OptVariableSet
{
  public:

    /**
     * Generic data struct for storing information about variables
     */
    struct { int type;
      int bus_id;
      int bus1_id;
      int bus2_id;
      double dmin;
      double dmax;
      int imin;
      int imax;
    } OptVar;

    /**
     * Constructor
     * @param network network from which variable set is derived
     */
    OptVariableSet(NetworkPtr network)
    {
      p_hasData = false;
    }

    /**
     * Destructor
     */
    ~OptVariableSet(void)
    {
      if (p_hasData) {
        GA_Deregister_type(p_gaType);
        GA_Destroy(p_gdata);
      }
    }

    /**
     * Initialize variable set by running over network and extracting all
     * variables from it
     */
    void init()
    {
      p_variables.clear();
      int numBus = p_network->numBuses();
      int numBranch = p_network->numBranches();
      int i, j, nvar, it;
      gridpack::component::OptimizationInterface *opt;
      // Loop over buses and branches and gather all variables
      for (i=0; i<numBus; i++) {
        if (p_network->getActiveBus(i)) {
          opt = dynamic_cast<gridpack::component::OptimizationInterface*>(
              p_network->getBus(i).get());
          nvar = opt->numOptVariables();
          OptVar var;
          for (j=0; j<nvar; j++) {
            it = opt->optVariableType(j);
            var.type = it;
            var.bus_id = p_network->getBus(i)->getOriginalBusIndex();
            var.bus1_id = -1;
            var.bus2_id = -1;
            if (it == OPT_INT) {
              opt->optVariableBounds(j,&var.imin,&var.imax);
              var.dmin = 0.0;
              var.dmax = 0.0;
            } else if (it == OPT_DBL) {
              opt->optVariableBounds(j,&var.dmin,&var.dmax);
              var.imin = 0;
              var.imax = 0;
            } else if (it == OPT_BOOL) {
              var.dmin = 0.0;
              var.dmax = 0.0;
              var.imin = 0;
              var.imax = 0;
            } else {
              // TODO: some kind of error
            }
            p_variables.push_back(var);
          }
        }
      }
      for (i=0; i<numBranch; i++) {
        if (p_network->getActiveBranch(i)) {
          opt = dynamic_cast<gridpack::component::OptimizationInterface*>(
              p_network->getBranch(i).get());
          nvar = opt->numOptVariables();
          OptVar var;
          for (j=0; j<nvar; j++) {
            it = opt->optVariableType(j);
            var.type = it;
            var.bus_id = 0;
            var.bus1_id = p_network->getBranch(i)->getOriginalBus1Index();
            var.bus2_id = p_network->getBranch(i)->getOriginalBus2Index();
            if (it == OPT_INT) {
              opt->optVariableBounds(j,&var.imin,&var.imax);
              var.dmin = 0.0;
              var.dmax = 0.0;
            } else if (it == OPT_DBL) {
              opt->optVariableBounds(j,&var.dmin,&var.dmax);
              var.imin = 0;
              var.imax = 0;
            } else if (it == OPT_BOOL) {
              var.dmin = 0.0;
              var.dmax = 0.0;
              var.imin = 0;
              var.imax = 0;
            } else {
              // TODO: some kind of error
            }
            p_variables.push_back(var);
          }
        }
      }
      // All variables from active buses and branches are now stored in
      // p_variables vector. Store data in global array
      if (p_hasData) {
        GA_Deregister_type(p_gaType);
        GA_Destroy(p_gdata);
      }
      p_gaType = NGA_Register_type(sizeof(OptVar));
      int nsize = p_variables.size();
      // Find total number of variables
      int nprocs = GA_Nnodes();
      int sizes[nprocs];
      int offsets[nprocs];
      int me = GA_Nodeid();
      for (i=0; i<nprocs; i++) {
        sizes[i] = 0;
      }
      sizes[me] = nsize;
      char cplus[2];
      strcpy(cplusm,"+");
      GA_Igop(sizes,nprocs,cplus);
      // Create global array
      p_totalVar = 0;
      for (i=0; i<nprocs; i++) {
        p_totalVar += sizes[i];
      }
      offsets[0] = 0;
      for (i=1; i<nprocs; i++) {
        offsets[i] = offsets[i-1]+sizes[i-1];
      }
      int one = 1;
      p_gdata = GA_Create_handle();
      NGA_Set_data(p_gdata,one,&p_totalVar,p_gaType);
      NGA_Set_irreg_distr(p_gdata,offsets,&nprocs);
      NGA_Allocate(p_gdata);
      // Allocate local buffer for variables
      OptVar *buf;
      buf = new OptVar[nsize];
      int lo = offsets[me];
      int hi = lo + nsize - 1;
      for (i=0; i<nsize; i++) {
        buf[i] = p_variables[i];
      }
      NGA_Put(p_gdata,&lo,&hi,buf,one);
      delete [] buf;
      GA_Sync();
    }

    /**
     * Return a list of all variables in the system from the OptVariableSet
     * object.
     * @return list of all variables
     */
    std::vector<OptVar> getFullVarList()
    {
      OptVar *buf;
      buf = new OptVar[p_totalVar];
      int lo = 0;
      int hi = p_totalVar-1;
      int one = 1;
      NGA_Get(p_gdata,&lo,&hi,buf,one);
      int i;
      std::vector<OptVar> ret;
      for (i=0; i<p_totalVar; i++) {
        ret.push_back(buf[i]);
      }
      delete [] buf;
      return ret;
    }

  private:
    NetworkPtr p_network;

    std::vector<OptVar> p_variables;

    int p_gaType;

    int p_gdata;

    bool p_hasData;

    int p_totalVar
};

}    // optimization
}    // gridpack
#endif
