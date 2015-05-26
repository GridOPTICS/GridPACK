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
//#include <gridpack/parallel/distributed.hpp>
#include <ga.h>


// Base optimization class that contains functions that are generic to all
// applications.

namespace gridpack{
namespace optimization{

template <class _network>
class Optimizer
//:
//  public gridpack::parallel::Communicator
{
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
    
    double *uc_iniLevel;
    double *uc_minUpTime ;
    double *uc_minDownTime ;
    double *uc_minPower ;
    double *uc_maxPower ;
    double *uc_costConst ;
    double *uc_costLinear ;
    double *uc_costQuad ;
    double *uc_rampUp ;
    double *uc_rampDown ;
    double *uc_startUp ;
    double *uc_initPeriod ;
    double *uc_startCap ;
    double *uc_shutCap ;
    double *uc_opMaxGen ;
    int totalGen;
    /**
     * Default Constructor
     * @param network - network the optimizer works on
     */
    Optimizer(NetworkPtr network)
      : p_network(network)
    { 
      p_nBuses = p_network->numBuses();
      totalGen = 0;
    }

//    Optimizer(const gridpack::parallel::Communicator& comm, NetworkPtr network)
//      : p_network(network),gridpack::parallel::Communicator()
//    { 
//      p_nBuses = p_network->numBuses();
//      totalGen = 0;
//    }

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
      * solution
      */
      
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

    /**
      * Get unit commitment parameters
      */  
    void getUCparam(void) 
    {
      int ngen;
      int gen_id;
//      gridpack::component::DataCollection *data;
      double rval;
// create global array to hold data for all generators in the network
      int grp = p_network->communicator().getGroup();
      int nprocs = GA_Pgroup_nnodes(grp);
      int me = GA_Pgroup_nodeid(grp);
      for (int i=0; i<p_nBuses; i++){
        if(p_network->getActiveBus(i)) {
          ngen = p_network->getBus(i)->numGen;
          totalGen += ngen;
        }
      }
      int loc_totalGen = totalGen;
      GA_Igop(&totalGen,1,"+");

      int genArr[nprocs];

      for (int p=0; p<nprocs; p++) {
        genArr[p] = 0;
      }
      genArr[me] = loc_totalGen;
      GA_Pgroup_igop(grp,genArr, nprocs, "+");
      int offset[nprocs];
      offset[0] = 0;
      for (int p=1; p<nprocs; p++) {
        offset[p]= offset[p-1] + genArr[p-1];
      }

      uc_iniLevel = new double[totalGen] ();
      uc_minUpTime = new double[totalGen] ();
      uc_minDownTime = new double[totalGen] ();
      uc_minPower = new double[totalGen] ();
      uc_maxPower = new double[totalGen] ();
      uc_costConst = new double[totalGen] ();
      uc_costLinear = new double[totalGen] ();
      uc_costQuad = new double[totalGen] ();
      uc_rampUp = new double[totalGen] ();
      uc_rampDown = new double[totalGen] ();
      uc_startUp = new double[totalGen] ();
      uc_initPeriod = new double[totalGen] ();
      uc_shutCap = new double[totalGen] ();
      uc_opMaxGen = new double[totalGen] ();
      uc_startCap = new double[totalGen] ();
      int index = offset[me];
      for (int i=0; i<p_nBuses; i++){
        if(p_network->getActiveBus(i)) {
          ngen = p_network->getBus(i)->numGen;
          if(ngen > 0) {
            for (int j=0; j<ngen; j++){
              uc_iniLevel[index] = p_network->getBus(i)->p_iniLevel[j];
              uc_minUpTime[index] = p_network->getBus(i)->p_minUpTime[j];
              uc_minDownTime[index] = p_network->getBus(i)->p_minDownTime[j];
              uc_minPower[index] = p_network->getBus(i)->p_minPower[j];
              uc_maxPower[index] = p_network->getBus(i)->p_maxPower[j];
              uc_costConst[index] = p_network->getBus(i)->p_costConst[j];
              uc_costLinear[index] = p_network->getBus(i)->p_costLinear[j];
              uc_costQuad[index] = p_network->getBus(i)->p_costQuad[j];
              uc_rampUp[index] = p_network->getBus(i)->p_rampUp[j];
              uc_rampDown[index] = p_network->getBus(i)->p_rampDown[j];
              uc_startUp[index] = p_network->getBus(i)->p_startUp[j];
              uc_shutCap[index] = p_network->getBus(i)->p_shutCap[j];
              uc_opMaxGen[index] = p_network->getBus(i)->p_opMaxGen[j];
              uc_initPeriod[index] = p_network->getBus(i)->p_initPeriod[j];
              uc_startCap[index] = p_network->getBus(i)->p_startCap[j];
              index++;
            }
          }
        }
      }
      GA_Pgroup_dgop(grp,uc_iniLevel, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_minUpTime, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_minDownTime, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_minPower, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_maxPower, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_costConst, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_costQuad, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_rampUp, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_rampDown, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_startUp, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_initPeriod, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_startCap, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_shutCap, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_costLinear, totalGen, "+");
      GA_Pgroup_dgop(grp,uc_opMaxGen, totalGen, "+");
    }
//
// print bus generator solution
    void solution(void) 
    {
      int numBus = p_network->numBuses();
      for (int i=0; i<numBus; i++) {
        if(p_network->getActiveBus(i)) {
          p_network->getBus(i)->solution();
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
