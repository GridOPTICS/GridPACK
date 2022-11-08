/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   gast.cpp
 * @author Shrirang Abhyankar
 * Added:  November 7, 2022
 * 
 * @brief  
 * Gas turbine governor model
 * 
 */

#include <vector>
#include <iostream>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_governor_model.hpp"
#include "gast.hpp"

#define TS_THRESHOLD 4

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::GastModel::GastModel(void)
{
  Pmech = 1.0;
  Pref = 1.0;
  delta_w = 0.0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::GastModel::~GastModel(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to
 * GastModel
 */
void gridpack::dynamic_simulation::GastModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(GOVERNOR_R, &R, idx)) R = 0.05; 
  if (!data->getValue(GOVERNOR_T1, &T1, idx)) T1 = 0.5;
  if (!data->getValue(GOVERNOR_T2, &T2, idx)) T2 = 3.0; 
  if (!data->getValue(GOVERNOR_T3, &T3, idx)) T3 = 10.0;
  if (!data->getValue(GOVERNOR_AT, &AT, idx)) AT = 0.0;
  if (!data->getValue(GOVERNOR_KT, &KT, idx))  KT = 0.0;
  if (!data->getValue(GOVERNOR_VMAX, &VMAX, idx)) VMAX = 1.0; 
  if (!data->getValue(GOVERNOR_VMIN, &VMIN, idx)) VMIN = 0.0; 
  if (!data->getValue(GOVERNOR_DT, &Dt, idx)) Dt = 0.0;

  fuel_valve_block.setparams(1.0,T1,VMIN,VMAX,-10000,10000);
  fuel_flow_block.setparams(1.0,T2);
  exh_temp_block.setparams(1.0,T3);
}

/**
 * Initialize governor model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::GastModel::init(double mag, double ang, double ts)
{
  //  Work backwards from Pmech
  double fuel_flow_block_out = Pmech + Dt*delta_w;
  
  double fuel_valve_block_out;

  // Initialize fuel valve block
  fuel_valve_block_out = fuel_flow_block.init_given_y(fuel_flow_block_out);

  // Initialize exh_temp_block
  exh_temp_block_out = exh_temp_block.init_given_u(fuel_flow_block_out);

  // Initialize fuel_valve_block
  double fuel_valve_block_in = fuel_valve_block.init_given_y(fuel_valve_block_out);
  double droop = delta_w/R;

  double lv_gate_in1 = AT + KT*(AT - exh_temp_block_out);

  Pref = fuel_valve_block_in + delta_w*R; // Here the assumption is lv_gate_in1 > fuel_valve_block_in.
  // If this condition is not satisfied then Pref value is indeterminate.
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GastModel::predictor(double t_inc, bool flag)
{
  double droop = delta_w/R;

  double lv_gate_in1 = AT + KT*(AT - exh_temp_block_out);

  double fuel_valve_block_in = std::min(lv_gate_in1,Pref - droop);

  // Fuel valve block
  double fuel_valve_block_out = fuel_valve_block.getoutput(fuel_valve_block_in,t_inc,PREDICTOR);

  // Fuel flow block
  double fuel_flow_block_out = fuel_flow_block.getoutput(fuel_valve_block_out,t_inc,PREDICTOR);

  Pmech = fuel_flow_block_out - Dt*delta_w;

  // Exhaust temperature block
  exh_temp_block_out = exh_temp_block.getoutput(fuel_flow_block_out,t_inc,PREDICTOR);
  
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GastModel::corrector(double t_inc, bool flag)
{
  double droop = delta_w/R;

  double lv_gate_in1 = AT + KT*(AT - exh_temp_block_out);

  double fuel_valve_block_in = std::min(lv_gate_in1,Pref - droop);

  // Fuel valve block
  double fuel_valve_block_out = fuel_valve_block.getoutput(fuel_valve_block_in,t_inc,CORRECTOR);

  // Fuel flow block
  double fuel_flow_block_out = fuel_flow_block.getoutput(fuel_valve_block_out,t_inc,CORRECTOR);

  Pmech = fuel_flow_block_out - Dt*delta_w;

  // Exhaust temperature block
  exh_temp_block_out = exh_temp_block.getoutput(fuel_flow_block_out,t_inc,CORRECTOR);

}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power
 */
void gridpack::dynamic_simulation::GastModel::setMechanicalPower(double pmech)
{
  Pmech = pmech; 
}

/**
 * Set the rotor speed deviation inside the governor
 * @param delta_o value of the rotor speed deviation
 */
void gridpack::dynamic_simulation::GastModel::setRotorSpeedDeviation(double delta_w)
{
  delta_w = delta_w;
}

/** 
 * Get the value of the mechanical power
 * @return value of mechanical power
 */
double gridpack::dynamic_simulation::GastModel::getMechanicalPower()
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
void gridpack::dynamic_simulation::GastModel::setExtBusNum(int ExtBusNum)
{
	p_bus_id = ExtBusNum;
}
*/	

/** 
 * Set the governor generator id
 */
 /**
void gridpack::dynamic_simulation::GastModel::setExtGenId(std::string ExtGenId)
{
	p_ckt = ExtGenId;
}
*/
