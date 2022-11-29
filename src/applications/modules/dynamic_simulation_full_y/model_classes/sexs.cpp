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

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::SexsModel::SexsModel(void)
{
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

  // Set parameters for the first block
  leadlagblock.setparams(TA,TB);

  zero_TE = false;
  if(fabs(TE) < 1e-6) {
    zero_TE = true;
  }
  
  // Set parameters for the second block
  if(!zero_TE) {
    filterblock.setparams(K,TE,EMIN,EMAX,-1000.0,1000);
  } else {
    gainblock.setparams(K,EMIN,EMAX);
  }

  Vs = 0.0; 
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::SexsModel::init(double Vm, double ang, double ts)
{
  double y1,u1;

  // For initialization, we are given the output for the model Efd.
  // To initialize the model blocks, we need to go backwards starting
  // from initializing second block and then the first one and then
  // calculating the model input Vref
  
  // Initialize second block
  if(!zero_TE) {
    y1 = filterblock.init_given_y(Efd);
  } else {
    y1 = std::min(EMAX,std::max(EMIN,Efd/K));
  }

  // Initialize first block
  u1 = leadlagblock.init_given_y(y1); 

  // Note: Vm is same as Ec
  Vref = Vm + u1 - Vs;   // Voltage reference initial value

}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::SexsModel::predictor(double t_inc, bool flag)
{
  double u1,y1;

  u1 = Vref + Vs - Ec;

  // Calculate first block output, last input flag = true
  // tells the block to do the state update (predictor update)
  y1 = leadlagblock.getoutput(u1,t_inc,PREDICTOR,true);

  // Calculate second block output, last input flag = true
  // tells the block to do the state update (predictor update)
  if(!zero_TE) {
    Efd = filterblock.getoutput(y1,t_inc,PREDICTOR,true);
  } else {
    Efd = gainblock.getoutput(y1);
  }
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::SexsModel::corrector(double t_inc, bool flag)
{
  double u1,y1;

  u1 = Vref + Vs - Ec;

  // Calculate first block output, last input flag = true
  // tells the block to do the state update (corrector update)
  y1 = leadlagblock.getoutput(u1,t_inc,CORRECTOR,true);

  // Calculate second block output, last input flag = true
  // tells the block to do the state update (corrector update)
  if(!zero_TE) {
    Efd = filterblock.getoutput(y1,t_inc,PREDICTOR,true);
  } else {
    Efd = gainblock.getoutput(y1);
  }
}

/**
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void gridpack::dynamic_simulation::SexsModel::setFieldVoltage(double fldv)
{
  // This is the initial value of Efd using during initialization
  Efd = fldv;
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double gridpack::dynamic_simulation::SexsModel::getFieldVoltage()
{
  double u1,y1;

  u1 = Vref + Vs - Ec;

  // Calculate first block output
  y1 = leadlagblock.getoutput(u1);

  // Calculate second block output
  if(!zero_TE) {
    Efd = filterblock.getoutput(y1);
  } else {
    Efd = gainblock.getoutput(y1);
  }

  return Efd;
}

/** 
 * Set the value of terminal voltage
 * 
 */
void gridpack::dynamic_simulation::SexsModel::setVterminal(double Vm)
{
  Ec = Vm;
}

void gridpack::dynamic_simulation::SexsModel::setVstab(double Vstab)
{
  Vs = Vstab;
}

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


