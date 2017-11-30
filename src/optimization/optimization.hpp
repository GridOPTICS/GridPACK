/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   optimization.hpp
 * @author 
 * @date   2015-11-20 11:52:32 d3g096
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
#include "gridpack/component/optimization_ifc.hpp"
#include <boost/smart_ptr/shared_ptr.hpp>

#include <ga.h>


// Base optimization class that contains functions that are generic to all
// applications.

namespace gridpack{
namespace optimization{

template <class _network>
class NetworkOptimizer
  : public gridpack::optimization::VariableVisitor
//:
//  public gridpack::parallel::Communicator
{
  public:
    typedef _network NetworkType;
    typedef boost::shared_ptr<NetworkType> NetworkPtr;
    typedef boost::shared_ptr<gridpack::optimization::Variable> VarPtr;
    typedef boost::shared_ptr<gridpack::optimization::Expression> ExpPtr;
    typedef boost::shared_ptr<gridpack::optimization::Constraint> ConstPtr;
    ExpPtr objFunc;
    std::vector<ConstPtr> locConstraint;
    

    int numUnits;
    std::vector<int> minUpTime;
    std::vector<int> minDownTime;
    std::vector<double> minPower;
    std::vector<double> demand;
    std::vector<double> reserve;
    std::vector<double> maxPower;
    std::vector<double> costConst;
    std::vector<double> costLinear;
    std::vector<double> costQuad;
    
    double *uc_iniLevel;
    double *uc_minUpTime ;
    double *uc_minDownTime ;
    double *uc_minPower ;
    double *uc_demand ;
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
    int *busID;
    int totalGen;
    /**
     * Default Constructor
     * @param network - network the optimizer works on
     */
    NetworkOptimizer(NetworkPtr network)
      : p_network(network)
    { 
      p_nBuses = p_network->numBuses();
      totalGen = 0;
    }


    /**
     * Destructor
     */
    ~NetworkOptimizer(void)
    {}

    /**
     * Return a vector of optimization variables associated witht this interface
     * @return list of variables
     */
    std::vector<boost::shared_ptr<gridpack::optimization::Variable> >
      getVariables()
    {
       std::vector<VarPtr> ret;
// to avoid accumulating
       VarPtr vptr;
       vptr->clear(); 
       ret.clear();
       double rval;
      int grp = p_network->communicator().getGroup();
      int nprocs = GA_Pgroup_nnodes(grp);
      int me = GA_Pgroup_nodeid(grp);
      int loc_totalGen = totalGen;

      std::vector<int> genArr(nprocs);

      for (int p=0; p<nprocs; p++) {
        genArr[p] = 0;
      }
      genArr[me] = p_numUnits*p_numHorizons*5;
      GA_Pgroup_igop(grp,&genArr[0], nprocs, "+");
      std::vector<int> offset(nprocs);
      offset[0] = 0;
      for (int p=1; p<nprocs; p++) {
        offset[p]= offset[p-1] + genArr[p-1];
      }
// 
// power produced
       for (int p = 0; p < p_numHorizons; p++) {
       int inc = offset[me];
         for (int i = 0; i < p_numUnits; i++) {
           rval = uc_maxPower[i];
           VarPtr vptr (new RealVariable(0.0, 0.0, rval));
           vptr->name(boost::str(boost::format("p_u_%d_t_%d")  %inc %p));

           ret.push_back(vptr);
           inc++;
         }
       }

// reserve required
       for (int p = 0; p < p_numHorizons; p++) {
       int inc = offset[me]; 
         for (int i = 0; i < p_numUnits; i++) {
           rval = uc_maxPower[i];
           VarPtr vptr (new RealVariable(0.0, 0.0, rval));
           vptr->name(boost::str(boost::format("r_u_%d_t_%d") %inc %p));
           ret.push_back(vptr);

           inc++;
         }
       }
// onoff 
       for (int p = 0; p < p_numHorizons; p++) {
       int inc = offset[me]; 
         for (int i = 0; i < p_numUnits; i++) {
           VarPtr vptr (new IntegerVariable(0,0,1));
           vptr->name(boost::str(boost::format("on_u_%d_t_%d") %inc %p));
           ret.push_back(vptr);

           inc++;
         }
       }
// start up 
       for (int p = 0; p < p_numHorizons; p++) {
       int inc = offset[me]; 
         for (int i = 0; i < p_numUnits; i++) {
           VarPtr vptr (new IntegerVariable(0,0,1));
           vptr->name(boost::str(boost::format("up_u_%d_t_%d") %inc %p));
           ret.push_back(vptr);

           inc++;
         }
       }
// shut down
       for (int p = 0; p < p_numHorizons; p++) {
       int inc = offset[me]; 
         for (int i = 0; i < p_numUnits; i++) {
           VarPtr vptr (new IntegerVariable(0,0,1));
           vptr->name(boost::str(boost::format("dn_u_%d_t_%d") %inc %p));
           ret.push_back(vptr);

           inc++;
         }
       }
       return ret;
    }

    /**
     * Return contribution from bus to a global constraint
     * @param tag string that can be parsed by bus to determine which constraint
     * contribution is being requested
     * @return contribution to global constraint. If no contribution, return
     * null pointer
     */
//    boost::shared_ptr<gridpack::optimization::Expression>
    std::vector<ExpPtr>
      getGlobalConstraint(const char* tag)
    {
      std::vector<VarPtr> vlist;
      vlist.clear();
      vlist = this->getVariables();
      int nVar = p_numHorizons*p_numUnits;
      int nVarP = nVar;
      VarPtr powerProduced;
      VarPtr powerReserved;
//  Initial state, treat as constraint
      int powerCnt;
      int reserveCnt;
// get the list of variables; clear() is required, otherwise, variable number is wrong
      std::vector<ExpPtr> ret;
      ret.clear();

// energy balance
      for (int p = 1; p < p_numHorizons; p++) {
        ExpPtr exprgP;
        exprgP.reset();
        ExpPtr exprgR;
        exprgR.reset();

        for (int i = 0; i < p_numUnits; i++) {
           powerCnt = p*p_numUnits + i;
           powerProduced = vlist[powerCnt];
           reserveCnt = nVarP + powerCnt;
           powerReserved = vlist[reserveCnt];
           if(!exprgP) {
             exprgP = 1*powerProduced;
           }else{
             exprgP = exprgP + 1*powerProduced;
           }

           if(!exprgR) {
             exprgR = 1*powerReserved;
           }else{
             exprgR = exprgR + 1*powerReserved;
           }
        }
        ret.push_back(exprgP);

        ret.push_back(exprgR);

      }
      return ret;
    }

    /**
     * Return a list of local constraints from component
     * @return list of constraints
     */
//    std::vector<boost::shared_ptr<gridpack::optimization::Constraint> >
    std::vector<ConstPtr>
      getLocalConstraints()
    {
      ConstPtr con;
      std::vector<VarPtr> vlist;
      vlist.clear();
      vlist = this->getVariables();
      ExpPtr expr;
      ExpPtr expr2;
      ExpPtr upDnIndicator;
      int nVar = p_numHorizons*p_numUnits;
      int nVarP = nVar;
      int nVarR = nVarP + nVar;
      int nIntVcntOnoff = nVarR + nVar;
      int nIntVcntStartup = nIntVcntOnoff + nVar;
      int nIntVcntShutdown = nIntVcntStartup + nVar;
      int powerCntm1 = 0;
      int onOffCntm1 = 0; 
      int onOffCnttmp = 0;
      int upDnPeriod;
      VarPtr onOff;
      VarPtr onOfftmp;
      VarPtr start_Up;
      VarPtr shutDown;
      VarPtr powerProduced;
      VarPtr powerReserved;
//  Initial state, treat as constraint
      int onOffCnt;
      int start_UpCnt;
      int shutDownCnt;
      int powerCnt;
      int reserveCnt;
// get the list of variables; clear() is required, otherwise, variable number is wrong
      std::vector<ConstPtr> ret;
      ret.clear();
        onOffCnt = nVarR + 1;
      for(int i=0; i< p_numUnits; i++) {
        onOffCnt = nVarR + i;
        con = vlist[onOffCnt] == 1;
        ret.push_back(con);

        start_UpCnt = nIntVcntOnoff + i;
        con = vlist[start_UpCnt] == 0;
        ret.push_back(con);

        shutDownCnt = nIntVcntStartup + i;
        con = vlist[shutDownCnt] == 0;
        ret.push_back(con);

        powerCnt = i;
        con = vlist[powerCnt] == uc_iniLevel[i];
        ret.push_back(con);

      }
      for (int p = 1; p < p_numHorizons; p++) {
        ExpPtr exprgP;
        exprgP.reset();
        ExpPtr exprgR;
        exprgR.reset();
        for (int i = 0; i < p_numUnits; i++) {
           powerCnt = p*p_numUnits + i;
           reserveCnt = nVarP + powerCnt;
           onOffCnt = nVarR + powerCnt;
           start_UpCnt = nIntVcntOnoff + powerCnt;
           shutDownCnt = nIntVcntStartup + powerCnt;
         
           onOff = vlist[onOffCnt];
           start_Up = vlist[start_UpCnt];
           shutDown = vlist[shutDownCnt];
           powerProduced = vlist[powerCnt];
           powerReserved = vlist[reserveCnt];

           expr = powerProduced - 10000*onOff;
           con = expr <= 0;
           ret.push_back(con);

           expr = powerProduced - uc_minPower[i]*onOff;
           con = expr >= 0;
           ret.push_back(con);

// ramp up constraint
           powerCntm1 = (p-1)*p_numUnits + i;
           expr = powerProduced + powerReserved - vlist[powerCntm1];
           con = expr <= uc_rampUp[i];
           ret.push_back(con);

// ramp down constraint
           expr = vlist[powerCntm1]-powerProduced;
           con = expr <= uc_rampDown[i];
           ret.push_back(con);

// minium up and down time
// on at horizon p
           onOffCntm1 = nVarR + (p-1)*p_numUnits + i;;
           upDnIndicator = onOff - vlist[onOffCntm1];
           if(p == 1) {
            int initP = (int)(uc_initPeriod[i]);
//
            upDnIndicator = vlist[onOffCntm1] - onOff;
            upDnPeriod = std::min(p_numHorizons, (p+(int)(uc_minUpTime[i]+0.5)-initP));
            for (int j = p; j < upDnPeriod; j++) {
              onOffCnttmp = nVarR + j*p_numUnits + i;;
              onOfftmp = vlist[onOffCnttmp];
              con = upDnIndicator - 10000*onOfftmp <= 0;
              ret.push_back(con);
            }
           } else{
             upDnPeriod = std::min(p_numHorizons, (p+(int)(uc_minUpTime[i]+0.5)));
            for (int j = p; j < upDnPeriod; j++) {
              onOffCnttmp = nVarR + j*p_numUnits + i;;
              onOfftmp = vlist[onOffCnttmp];
              con = upDnIndicator - 10000*onOfftmp <= 0;
              ret.push_back(con);
            }
           }

// start up, previous off
           upDnIndicator = onOff - vlist[onOffCntm1];
           con = upDnIndicator - 10000*start_Up <= 0;
           ret.push_back(con);

// off at horizon p
           upDnIndicator = - 1*vlist[onOffCntm1] + 1*onOff;
           upDnPeriod = std::min(p_numHorizons, (p+(int)(uc_minDownTime[i]+0.5)));
           for (int j = p; j < upDnPeriod; j++) {
             onOffCnttmp = nVarR + j*p_numUnits + i;;
             onOfftmp = vlist[onOffCnttmp];
             con = upDnIndicator - 10000*onOfftmp <= -1;
             ret.push_back(con);
           }
// shut down, previous on
           upDnIndicator = vlist[onOffCntm1] - onOff;
           con = upDnIndicator - 10000*shutDown <= 0;
           ret.push_back(con);

// generation limits
// startup at horizon p
           expr = powerProduced + powerReserved
           - uc_maxPower[i]*onOff + uc_minPower[i]*onOff
               + uc_maxPower[i]*start_Up - uc_startCap[i]*start_Up;
           con = expr <= uc_minPower[i];
           ret.push_back(con);

// shutdown at horizon p 
           onOffCntm1 = nVarP + powerCntm1;
           onOffCnttmp = nVarR + powerCntm1;
           expr = vlist[powerCntm1] + vlist[onOffCntm1] -
            uc_maxPower[i]*vlist[onOffCnttmp] + uc_minPower[i]*vlist[onOffCnttmp]
                +uc_maxPower[i]*shutDown - uc_shutCap[i]*shutDown;
           con = expr  <= uc_minPower[i];
           ret.push_back(con);
        }

      }
      return ret;
    }


    /**
     * Return contribution to objective function
     * @return expression representing contribution to objective function. If no
     * contribution, return null pointer
     */
//    boost::shared_ptr<gridpack::optimization::Expression>
    ExpPtr
      getObjectiveFunction()
    {
//  boost::shared_ptr<gridpack::optimization::Expression> ret;
      int nVar = p_numHorizons*p_numUnits;
      int nRealVcnt = 0;
      int nIntVcnt = 2*nVar;
      std::vector<VarPtr> vlist;
      vlist.clear();
      vlist = this->getVariables();
      ExpPtr obj;
      VarPtr onOff;
      VarPtr powerProduced;
      VarPtr start_Up;
      for (int p = 1; p < p_numHorizons; p++) {
        for (int i = 0; i < p_numUnits; i++) {
           nRealVcnt = p*p_numUnits + i;
           nIntVcnt = nRealVcnt + 2*nVar;
           powerProduced = vlist[nRealVcnt];
           onOff = vlist[nIntVcnt];
           start_Up = vlist[nIntVcnt+nVar];
           
         if(!obj) {
           obj =  uc_costConst[i]*onOff
              + uc_startUp[i]*start_Up
              + uc_costLinear[i]*powerProduced
              + 2.0*uc_costQuad[i]*((powerProduced)^2);
         }else
         {
           obj = obj + uc_costConst[i]*onOff
              + uc_startUp[i]*start_Up
              + uc_costLinear[i]*powerProduced
              + 2.0*uc_costQuad[i]*((powerProduced)^2);
         }
        }
      }
      return obj;
    }


    /**
     * sum over processes to get global objective function
     */
    double objectiveFunction(void)
    {
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
          data->getValue("GENERATOR_DEMAND", &rvar, idx);
          demand.push_back(rvar);
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
     * load bus data for expression test 
     */
    void loadBusData_exp(void) 
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
//      GA_Igop(&totalGen,1,"+");

      std::vector<int> genArr(nprocs);

      for (int p=0; p<nprocs; p++) {
        genArr[p] = 0;
      }
      genArr[me] = loc_totalGen;
      GA_Pgroup_igop(grp,&genArr[0], nprocs, "+");
      std::vector<int> offset(nprocs);
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
      busID = new int[totalGen] ();
//      int index = offset[me];
      int index = 0;
      int busid;
      for (int i=0; i<p_nBuses; i++){
        if(p_network->getActiveBus(i)) {
          ngen = p_network->getBus(i)->numGen;
          if(ngen > 0) {
            for (int j=0; j<ngen; j++){
              p_network->getBusData(i)->getValue("GENERATOR_BUSNUMBER",&busid,j);
              busID[index] = busid;        
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
#if 0
      GA_Pgroup_igop(grp,busID, totalGen, "+");
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
#endif
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
// Get number of periods and generators
    void getLoadsInfo(int numHorizons,double *demand_in, double *reserve_in)
    {
      int ngen;
      p_numUnits = 0;
      for (int i=0; i<p_nBuses; i++){
        if(p_network->getActiveBus(i)) {
          ngen = p_network->getBus(i)->numGen;
          p_numUnits += ngen;
        }
      }
      p_numHorizons = numHorizons;
//
      demand.clear();
      reserve.clear();
      for (int i=0; i<numHorizons; i++) {
        demand.push_back(demand_in[i]);
        reserve.push_back(reserve_in[i]);
      }
    }

  protected:

    NetworkPtr p_network;

  private:
    int  p_nBuses;
    int p_numUnits;
    int p_numHorizons;
};

}    // optimization
}    // gridpack
#endif
