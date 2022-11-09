/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   hygov.cpp
 * @author Shrirang Abhyankar
 * Added:  November 7, 2022
 * 
 * @brief  
 * Hydro turbine governor model
 * 
 */

#include <vector>
#include <iostream>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_governor_model.hpp"
#include "hygov.hpp"

#define TS_THRESHOLD 4

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::HygovModel::HygovModel(void)
{
  Pmech = 1.0;
  nref = 1.0;
  delta_w = 0.0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::HygovModel::~HygovModel(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to
 * HygovModel
 */
void gridpack::dynamic_simulation::HygovModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(GOVERNOR_R, &R, idx)) R = 0.05;
  if (!data->getValue(GOVERNOR_r, &r, idx)) R = 0.05; 
  if (!data->getValue(GOVERNOR_TR, &TR, idx)) TR = 0.5;
  if (!data->getValue(GOVERNOR_TF, &TF, idx)) TF = 3.0; 
  if (!data->getValue(GOVERNOR_TG, &TG, idx)) TG = 10.0;
  if (!data->getValue(GOVERNOR_VELM, &VELM, idx)) VELM = 1000.0;
  if (!data->getValue(GOVERNOR_GMAX, &GMAX, idx)) GMAX = 1.0; 
  if (!data->getValue(GOVERNOR_GMIN, &GMIN, idx)) GMIN = 0.0; 
  if (!data->getValue(GOVERNOR_TW, &TW, idx)) TW = 10.0;
  if (!data->getValue(GOVERNOR_AT, &AT, idx)) AT = 1.0;
  if (!data->getValue(GOVERNOR_DT, &Dt, idx)) Dt = 0.0;
  if (!data->getValue(GOVERNOR_QNL, &qNL, idx))  qNL = 0.0;

  filter_block.setparams(1.0,TF);
  gate_block.setparams(1/r,1/(r*TR),GMIN,GMAX,-10000,10000);
  opening_block.setparams(1.0,TG);
  turbine_flow_block.setparams(TW);
}

/**
 * Initialize governor model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::HygovModel::init(double mag, double ang, double ts)
{
  turbine_flow_block_out = Pmech/AT + qNL;

  double turbine_flow_block_in;

  turbine_flow_block_in = turbine_flow_block.init_given_y(turbine_flow_block_out);

  opening_block_out = turbine_flow_block_out;

  gate_block_out = opening_block.init_given_y(opening_block_out);

  double gate_block_in;
  gate_block_in = gate_block.init_given_y(gate_block_out);
  
  filter_block_in = filter_block.init_given_y(gate_block_in);

  nref = filter_block_in + R*gate_block_out + delta_w;
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::HygovModel::predictor(double t_inc, bool flag)
{
  filter_block_in = nref - (delta_w + R*gate_block_out);

  double filter_block_out = filter_block.getoutput(filter_block_in,t_inc,PREDICTOR,true);

  gate_block_out = gate_block.getoutput(filter_block_out,t_inc,PREDICTOR,true);

  opening_block_out = opening_block.getoutput(gate_block_out,t_inc,PREDICTOR,true);

  double turbine_flow_block_in,h;

  h = opening_block_out/turbine_flow_block_out;
  h = h*h;

  turbine_flow_block_in = 1 - h;

  turbine_flow_block_out = turbine_flow_block.getoutput(turbine_flow_block_in,t_inc,PREDICTOR,true);

  Pmech =  AT*(turbine_flow_block_out - qNL)*h - opening_block_out*Dt*delta_w;
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::HygovModel::corrector(double t_inc, bool flag)
{
  filter_block_in = nref - (delta_w + R*gate_block_out);

  double filter_block_out = filter_block.getoutput(filter_block_in,t_inc,CORRECTOR,true);

  gate_block_out = gate_block.getoutput(filter_block_out,t_inc,CORRECTOR,true);

  opening_block_out = opening_block.getoutput(gate_block_out,t_inc,CORRECTOR,true);

  double turbine_flow_block_in,h;

  h = opening_block_out/turbine_flow_block_out;
  h = h*h;

  turbine_flow_block_in = 1 - h;

  turbine_flow_block_out = turbine_flow_block.getoutput(turbine_flow_block_in,t_inc,CORRECTOR,true);

  Pmech =  AT*(turbine_flow_block_out - qNL)*h - opening_block_out*Dt*delta_w;

}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power
 */
void gridpack::dynamic_simulation::HygovModel::setMechanicalPower(double pmech)
{
  Pmech = pmech; 
}

/**
 * Set the rotor speed deviation inside the governor
 * @param delta_o value of the rotor speed deviation
 */
void gridpack::dynamic_simulation::HygovModel::setRotorSpeedDeviation(double delta_w)
{
  delta_w = delta_w;
}

/** 
 * Get the value of the mechanical power
 * @return value of mechanical power
 */
double gridpack::dynamic_simulation::HygovModel::getMechanicalPower()
{
  return Pmech; 
}

/** 
 * Get the value of the rotor speed deviation
 * 
 */

/** 
 * Set the governor generator bus number
 */
 /**
void gridpack::dynamic_simulation::HygovModel::setExtBusNum(int ExtBusNum)
{
	p_bus_id = ExtBusNum;
}
*/	

/** 
 * Set the governor generator id
 */
 /**
void gridpack::dynamic_simulation::HygovModel::setExtGenId(std::string ExtGenId)
{
	p_ckt = ExtGenId;
}
*/
