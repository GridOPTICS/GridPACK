/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   base_governor_model.cpp
 * @author Bruce Palmer
 * @Last modified:   May 18, 2015
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_governor_model.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::BaseGovernorModel::BaseGovernorModel(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::BaseGovernorModel::~BaseGovernorModel(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to
 * BaseGovernorModel
 */
void gridpack::dynamic_simulation::BaseGovernorModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
}

/**
 * Initialize governor model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::BaseGovernorModel::init(
    double mag, double ang, double ts)
{
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseGovernorModel::predictor(
    double t_inc, bool flag)
{
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseGovernorModel::corrector(
    double t_inc, bool flag)
{
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power
 */
void gridpack::dynamic_simulation::BaseGovernorModel::
setMechanicalPower(double pmech)
{
}

/**
 * Set the rotor speed deviation inside the governor
 * @param delta_o value of the rotor speed deviation
 */
void gridpack::dynamic_simulation::BaseGovernorModel::
setRotorSpeedDeviation(double delta_o)
{
}

/** 
 * Get the value of the mechanical power
 * @return value of mechanical power
 */
double gridpack::dynamic_simulation::BaseGovernorModel::
getMechanicalPower()
{
  return 0.0;
}

/** 
 * Get the value of the rotor speed deviation
 * 
 */
double gridpack::dynamic_simulation::BaseGovernorModel::
getRotorSpeedDeviation()
{
  return 0.0;
}
