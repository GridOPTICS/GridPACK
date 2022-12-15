/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   wtara1.cpp
 * @author Shuangshuang Jin 
 * @Created on:    Nov 15, 2022
 * 
 * @Last Updated: Dec 7, 2022
 * Shrirang Abhyankar
 * Added all the pieces
 *
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <string>
#include <cstdio>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "wtara1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Wtara1Model::Wtara1Model(void)
{
    Theta = 0.0;
    Pmech0 = 0.0;
    domega_t = 0.0;
    Taero = 0.0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Wtara1Model::~Wtara1Model(void)
{
}

/**
 * Load parameters from DataCollection object into mechanical model
 * @param data collection of mechanical parameters from input files
 * @param index of mechanical on bus
 * TODO: might want to move this functionality to
 * Wtara1Model
 */
void gridpack::dynamic_simulation::Wtara1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(WIND_AD_KA, &Ka, idx)) Ka = 0.007; // Ka
  if (!data->getValue(WIND_AD_THETA, &Theta0, idx)) Theta0 = 0.0; // theta
}

/**
 * Initialize  model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Wtara1Model::init(double mag, double ang, double ts)
{
  double Pmech;
  Pmech = Pmech0 - Ka * Theta * (Theta - Theta0);

  Taero = Pmech/(1 + domega_t);

}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtara1Model::predictor(double t_inc, bool flag)
{
  double Pmech;
  Pmech = Pmech0 - Ka * Theta * (Theta - Theta0);

  Taero = Pmech/(1 + domega_t);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtara1Model::corrector(double t_inc, bool flag)
{
  double Pmech;
  Pmech = Pmech0 - Ka * Theta * (Theta - Theta0);

  Taero = Pmech/(1 + domega_t);
}

/**
 * Set the value of initial mechanical power
 * @param: Pmech0 => Initial mechanical power
 */
void  gridpack::dynamic_simulation::Wtara1Model::setPmech(double Pmech0_in)
{
  Pmech0 = Pmech0_in;
}

/**
 * setTurbineSpeedDeviation - sets the turbine speed deviation
 * @param domega_turb : turbine speed deviation
 * From drive train model
 **/
void gridpack::dynamic_simulation::Wtara1Model::setTurbineSpeedDeviation(double domega_turb)
{
  domega_t = domega_turb;
}

/**
 * setTheta - sets pitch angle
 * @param Theta - pitch angle
 **/
void gridpack::dynamic_simulation::Wtara1Model::setTheta(double Theta_in)
{
  Theta = Theta_in;
}

/**
 * getTaero - returns the aero-dynamic torque
 * @output Taero - aero dynamic torque, output of aerodynamic model
 **/
double gridpack::dynamic_simulation::Wtara1Model::getTaero()
{
  return Taero;
}

/**
 * getTheta - Get initial value of pitch controller
 * @output theta0 - initial value of pitch controller
 **/
double gridpack::dynamic_simulation::Wtara1Model::getTheta() {
  return Theta0;
}
