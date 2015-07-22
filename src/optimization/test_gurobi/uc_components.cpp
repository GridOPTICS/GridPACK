/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   uc_components.cpp
 * @author 
 * @date   
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "uc_components.hpp"
#include "gridpack/parser/dictionary.hpp"

//#define LARGE_MATRIX

/**
 *  Simple constructor
 */
gridpack::unit_commitment::UCBus::UCBus(void)
{
  p_num_generator = 0;
  numGen = 0;
  p_bus_id = 0;
}

/**
 *  Simple destructor
 */
gridpack::unit_commitment::UCBus::~UCBus(void)
{
}


/**
 * Load values stored in DataCollection object into UCBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::unit_commitment::UCBus::load(const
         boost::shared_ptr<gridpack::component::DataCollection> &data)
{
   double rval;
   int gen_id;
   int ngen;
//   data->getValue("BUS_ID",&p_bus_id,i);
   if(data->getValue("GENERATOR_NUMBER", &ngen)) {
     numGen += ngen;
     for( int i=0; i< ngen; i++) {
       data->getValue("GENERATOR_ID",&gen_id,i);
       data->getValue("GENERATOR_INIT_LEVEL", &rval,i);
       p_iniLevel.push_back(rval);
       data->getValue("GENERATOR_MIN_GEN",&rval,i);
//printf ("blevel %f %f %d\n",rval,rval,ngen);
       p_minPower.push_back(rval);
       data->getValue("GENERATOR_MAX_GEN",&rval,i);
       p_maxPower.push_back(rval);
       data->getValue("GENERATOR_MIN_UP",&rval,i);
       p_minUpTime.push_back(rval);
       data->getValue("GENERATOR_MIN_DOWN",&rval,i);
       p_minDownTime.push_back(rval);
       data->getValue("GENERATOR_CONST_COST",&rval,i);
       p_costConst.push_back(rval);
       data->getValue("GENERATOR_LIN_COST",&rval,i);
       p_costLinear.push_back(rval);
       data->getValue("GENERATOR_CO_2_COST",&rval,i);
       p_costQuad.push_back(rval);
       data->getValue("GENERATOR_MAX_OP_GEN",&rval,i);
       p_opMaxGen.push_back(rval) ;  
       data->getValue("GENERATOR_RAMP_UP",&rval,i);
       p_rampUp.push_back(rval) ;
       data->getValue("GENERATOR_RAMP_DOWN",&rval,i);
       p_rampDown.push_back(rval) ;
       data->getValue("GENERATOR_START_UP",&rval,i);
       p_startUp.push_back(rval) ;
       data->getValue("GENERATOR_START_CAP",&rval,i);
       p_startCap.push_back(rval) ;
       data->getValue("GENERATOR_SHUT_CAP",&rval,i);
       p_shutCap.push_back(rval) ;
       data->getValue("GENERATOR_INIT_PRD",&rval,i);
       p_initPeriod.push_back(rval) ;
      }
   }
}
/**
  * objective cost
  */
double objectiveFunction(void) {
  double obj;
  obj = 2.0;
  return obj;
}
/**
  * uc solution
  */
bool gridpack::unit_commitment::UCBus::solution(void)
{
   double rval;
#if 0
   if(data->getValue("GENERATOR_NUMBER", &ngen)) {
     for( int i=0; i< ngen; i++) {
       data->getValue("POWER_PRODUCED", &rval,i);
       p_powerProduced.push_back(rval);
       data->getValue("POWER_RESERVED", &rval,i);
       p_powerReserved.push_back(rval);
    }
  }
#endif
  return true;
}

/**
 *  Simple constructor
 */
gridpack::unit_commitment::UCBranch::UCBranch(void)
{
}

/**
 *  Simple destructor
 */
gridpack::unit_commitment::UCBranch::~UCBranch(void)
{
}

/**
 * Load values stored in DataCollection object into UCBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::unit_commitment::UCBranch::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
}

