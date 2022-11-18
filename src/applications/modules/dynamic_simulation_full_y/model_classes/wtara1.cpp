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
#include "base_mechanical_model.hpp"
#include "wtara1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Wtara1Model::Wtara1Model(void)
{
    Theta = 0.0;
    Pmech = 0.0;
    wt = 0.0;
    waerot = 0.0;
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
    Pm0 = 1.0;
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtara1Model::predictor(double t_inc, bool flag)
{
    Pmech = Pm0 - Ka * Theta * (Theta - Theta0);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtara1Model::corrector(double t_inc, bool flag)
{
    Pmech = Pm0 - Ka * Theta * (Theta - Theta0);
}

/**
 * Set the Theta parameter inside the Aerodynamic model
 * @param Theta value of the pitch angle
 */
void  gridpack::dynamic_simulation::Wtara1Model::setPitchAngle(double Theta1)
{
    Theta = Theta1;
}

/**
 * Get the value of the mechanical power
 * @return value of mechanical power
 */
double  gridpack::dynamic_simulation::Wtara1Model::getMechanicalPower()
{
  return Pmech;
}

/**
 * Set the wt parameter inside the Aerodynamic model
 * @param wt value of the turbine speed
 */
void  gridpack::dynamic_simulation::Wtara1Model::setWt(double wt1)
{
    wt = wt1;
}

/**
 * Get the waerot parameter inside the Aerodynamic model
 * @return waerot value of the turbine speed
 */
double  gridpack::dynamic_simulation::Wtara1Model::getWaerot()
{
    waerot = Pmech / wt;
    return waerot;
}
