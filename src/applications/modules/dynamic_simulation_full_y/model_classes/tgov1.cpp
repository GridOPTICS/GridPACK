/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   tgov1.cpp
 * @author Shrirang Abhyankar
 * @Last modified:   November 29, 2022
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_governor_model.hpp"
#include "tgov1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Tgov1Model::Tgov1Model(void)
{
  Pmech = 1.0;
  Pref = 1.0;
  delta_w = 0.0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Tgov1Model::~Tgov1Model(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to
 * Tgov1Model
 */
void gridpack::dynamic_simulation::Tgov1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(GOVERNOR_R, &R, idx)) R = 0.05; 
  if (!data->getValue(GOVERNOR_T1, &T1, idx)) T1 = 0.5; 
  if (!data->getValue(GOVERNOR_VMAX, &Vmax, idx)) Vmax = 1.0; 
  if (!data->getValue(GOVERNOR_VMIN, &Vmin, idx)) Vmin = 0.0; 
  if (!data->getValue(GOVERNOR_T2, &T2, idx)) T2 = 3.0; 
  if (!data->getValue(GOVERNOR_T3, &T3, idx)) T3 = 10.0; 
  if (!data->getValue(GOVERNOR_DT, &Dt, idx)) Dt = 0.0;

  // Set up transfer function blocks
  leadlag_blk.setparams(T2,T3);
  delay_blk.setparams(1.0,T1,Vmin,Vmax,-1000.0,1000.0);
}

/**
 * Initialize governor model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Tgov1Model::init(double mag, double ang, double ts)
{
  // Initialize leadlag block
  delay_blk_out = leadlag_blk.init_given_y(Pmech+Dt*delta_w);

  double delay_blk_in;

  // Initialize delay block
  delay_blk_in = delay_blk.init_given_y(delay_blk_out);

  // Reference power signal
  Pref = R*delay_blk_in + delta_w;
}

void gridpack::dynamic_simulation::Tgov1Model::computeModel(double t_inc, IntegrationStage int_flag)
{
  // Delay block output and state update
  delay_blk_out = delay_blk.getoutput((Pref-delta_w)/R,t_inc,int_flag,true);

  // Leadlag block output and state update
  leadlag_blk_out = leadlag_blk.getoutput(delay_blk_out,t_inc,int_flag,true);

  // Output mechanical power
  Pmech = leadlag_blk_out - Dt*delta_w;
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Tgov1Model::predictor(double t_inc, bool flag)
{
  computeModel(t_inc,PREDICTOR);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Tgov1Model::corrector(double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR);
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power
 */
void gridpack::dynamic_simulation::Tgov1Model::setMechanicalPower(double Pmech0)
{
  Pmech = Pmech0; 
}

/**
 * Set the rotor speed deviation inside the governor
 * @param delta_o value of the rotor speed deviation
 */
void gridpack::dynamic_simulation::Tgov1Model::setRotorSpeedDeviation(double dw)
{
  delta_w = dw;
}

/** 
 * Get the value of the mechanical power
 * @return value of mechanical power
 */
double gridpack::dynamic_simulation::Tgov1Model::getMechanicalPower()
{
  return Pmech; 
}

/**
 * Set internal state parameter in governor
 * @param name character string corresponding to state variable
 * @param value new value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::Tgov1Model::setState(std::string name,
    double value)
{
  return false;
}

/**
 * Get internal state parameter in governor
 * @param name character string corresponding to state variable
 * @param value current value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::Tgov1Model::getState(std::string name,
    double *value)
{
  return false;
}
