/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   optimization.hpp
 * @author 
 * @date   
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _optimization_hpp_
#define _optimization_hpp_

#include <vector>
#include <macdecls.h>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"

// Base optimization class that contains functions that are generic to all
// applications.

namespace gridpack{
namespace optimization{

template <class _network>
class Optimizer {
  public:
    typedef _network NetworkType;
    typedef boost::shared_ptr<NetworkType> NetworkPtr;
    int numUnits;
    std::vector<int> minUpTime;
    std::vector<int> minDownTime;
    std::vector<double> minPower;
    std::vector<double> maxPower;
    std::vector<double> costConst;
    std::vector<double> costLinear;
    std::vector<double> costQuad;

    /**
     * Constructor
     * @param network - network the optimizer works on
     */
    Optimizer(NetworkPtr network)
      : p_network(network)
    { 
      p_nBuses = p_network->numBuses();
    }

    /**
     * Destructor
     */
    ~Optimizer(void)
    {}

    /**
     * sum over processes to get global objective function
     */
    double objectiveFunction(void)
    {
//      gridpack::utility::CoarseTimer *timer =
//        gridpack::utility::CoarseTimer::instance();
//      timer->configTimer(p_profile);
//      int t_setc = timer->createCategory("Factory:setComponents");
//      timer->start(t_setc);
      int numBus = p_network->numBuses();
      int i, j;
      int idx1, idx2;
      double sum;
      sum = 0;
//    loop over buses;
      for (i=0; i<numBus; i++) {
        if(p_network->getActiveBus(i)) {
          int branch_idx, bus1_idx, bus2_idx;
          sum += p_network->getBus(i)->objectiveFunction();
        }
      }
      GA_Dgop(&sum,1,"+");
      return sum;
//      timer->stop(t_setc);
//      timer->configTimer(true);
    }

    /**
     * load bus data 
     */
    void loadBusData(void) 
    {
      int ivar;
      double rvar;
      gridpack::component::DataCollection *data;
      for (int i=0; i<p_nBuses; i++) {
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(i).get());
        data->getValue("GENERATOR_NUMBERS",&numUnits);
// To sum to get total number of generators on the network
        for (int idx = 0; idx<numUnits; idx++) {
          data->getValue("GENERATOR_MIN_POWER", &rvar, idx);
          minPower.push_back(rvar);
        }
        for (int idx = 0; idx<numUnits; idx++) {
          data->getValue("GENERATOR_MAX_POWER", &rvar, idx);
          maxPower.push_back(rvar);
        }
        for (int idx = 0; idx<numUnits; idx++) {
          data->getValue("GENERATOR_MIN_UPTIME", &ivar, idx);
          minUpTime.push_back(ivar);
        }
        for (int idx = 0; idx<numUnits; idx++) {
          data->getValue("GENERATOR_MIN_DNTIME", &ivar, idx);
          minDownTime.push_back(ivar);
        }
        for (int idx = 0; idx<numUnits; idx++) {
          data->getValue("GENERATOR_CONST_COST", &rvar, idx);
          costConst.push_back(rvar);
        }
        for (int idx = 0; idx<numUnits; idx++) {
          data->getValue("GENERATOR_LINEAR_COST", &rvar, idx);
          costLinear.push_back(rvar);
        }
        for (int idx = 0; idx<numUnits; idx++) {
          data->getValue("GENERATOR_QUAD_COST", &rvar, idx);
          costQuad.push_back(rvar);
        }
      } 
    }

  protected:

    NetworkPtr p_network;

  private:
    int  p_nBuses;
};

}    // optimization
}    // gridpack
#endif
