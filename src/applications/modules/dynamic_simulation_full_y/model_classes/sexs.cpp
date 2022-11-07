/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   sexs.cpp
 * @author Shrirang Abhyankar
 * @Added:   Nov 6, 2022
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <cstdio>

// Yuan added below 2020-6-23
//#include <cstdio>
#include <cstring>
#include <string>
// Yuan added above 2020-6-23

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_exciter_model.hpp"
#include "sexs.hpp"

#define TS_THRESHOLD 4

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::SexsModel::SexsModel(void)
{
  OptionToModifyLimitsForInitialStateLimitViolation = true;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::SexsModel::~SexsModel(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * SexsModel
 */
void gridpack::dynamic_simulation::SexsModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(EXCITER_TA_OVER_TB, &TA_OVER_TB, idx)) TA_OVER_TB = 0.0; // TA_OVER_TB
  if (!data->getValue(EXCITER_TB, &TB, idx)) TB = 0.0; // TB
  if (!data->getValue(EXCITER_K, &K, idx))   K  = 0.0; // K
  if (!data->getValue(EXCITER_EMAX, &EMAX, idx)) EMAX = 0.0; // EMAX
  if (!data->getValue(EXCITER_EMIN, &EMIN, idx)) EMIN = 0.0; // EMIN
  if (!data->getValue(EXCITER_TE, &TE, idx)) TE = 0.0; // TE

  TA = TA_OVER_TB*TB;

  leadlagblock.setparams(TA,TB);
  filterblock.setparams(K,TE,EMIN,EMAX,-1000.0,1000);
}

/**
 *  * Saturation function
 *   * @ param x
 *    */
double gridpack::dynamic_simulation::SexsModel::Sat(double x)
{
	return 0;
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::SexsModel::init(double mag, double ang, double ts)
{
  filterblock.init(0.0,Efd); // Initialize second integrator
  double y1 = Efd/K;         // Output of first block

  leadlagblock.init(0.0,y1); // Initialize first integrator
  double u1 = y1/(1 - TA + TA_OVER_TB); // Input for first integrator

  Vstab = 0.0; // No stabilizer signal implemented
  
  Vref = mag + u1;   // Voltage reference

}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::SexsModel::predictor(double t_inc, bool flag)
{
  double u1,y1;

  u1 = Vref - Vterminal;
  
  y1 = leadlagblock.getoutput(u1,t_inc,PREDICTOR);

  Efd = filterblock.getoutput(y1,t_inc,PREDICTOR);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::SexsModel::corrector(double t_inc, bool flag)
{
  double u1,y1;

  u1 = Vref - Vterminal;
  
  y1 = leadlagblock.getoutput(u1,t_inc,CORRECTOR);

  Efd = filterblock.getoutput(y1,t_inc,CORRECTOR);
}

/**
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void gridpack::dynamic_simulation::SexsModel::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

/**
 * Set the field current parameter inside the exciter
 * @param fldc value of the field current
 */
void gridpack::dynamic_simulation::SexsModel::setFieldCurrent(double fldc)
{
  LadIfd = fldc;
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double gridpack::dynamic_simulation::SexsModel::getFieldVoltage()
{
  return Efd;
}

/** 
 * Get the value of the field current parameter
 * @return value of field current
 */
double gridpack::dynamic_simulation::SexsModel::getFieldCurrent()
{
  return 0.0;
}

/** 
 * Set the value of the Vterminal
 * @return value of field current
 */
void gridpack::dynamic_simulation::SexsModel::setVterminal(double mag)
{
  Vterminal = mag;
}

/** 
 * Set the value of the omega
 * @return value of field current
 */
void gridpack::dynamic_simulation::SexsModel::setOmega(double omega)
{
}

void gridpack::dynamic_simulation::SexsModel::setVstab(double vtmp)
{
  Vstab = vtmp;
}


// Yuan added below 2020-6-23
/** 
 * Set the exciter bus number
 * @return value of exciter bus number
 */
void gridpack::dynamic_simulation::SexsModel::setExtBusNum(int ExtBusNum)
{
	p_bus_id = ExtBusNum;
}	

/** 
 * Set the exciter generator id
 * @return value of generator id
 */
void gridpack::dynamic_simulation::SexsModel::setExtGenId(std::string ExtGenId)
{
	p_ckt = ExtGenId;
}	
// Yuan added above 2020-6-23

