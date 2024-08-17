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
 * Update parameters in DataCollection object with current values from
 * governor
 * @param data collection object for bus that hosts governor
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::BaseGovernorModel::updateData(
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

// Yuan added below 2020-6-23
/** 
 * Set the governor bus number
 */
void gridpack::dynamic_simulation::BaseGovernorModel::setExtBusNum(int ExtBusNum)
{
}	

/** 
 * Set the governor generator id
 */
void gridpack::dynamic_simulation::BaseGovernorModel::setExtGenId(std::string ExtGenId)
{
}	
// Yuan added above 2020-6-23

/**
 * Set internal state parameter in governor
 * @param name character string corresponding to state variable
 * @param value new value for state parameter
 * @return false if no variable corresponding to name is
 found
 */
bool gridpack::dynamic_simulation::BaseGovernorModel::setState(
    std::string name, double value)
{
  return false;
}

/**
 * Get internal state parameter in governor
 * @param name character string corresponding to state variable
 * @param value current value for state parameter
 * @return false if no variable corresponding to name is
 found
 */
bool gridpack::dynamic_simulation::BaseGovernorModel::getState(
    std::string name, double *value)
{
  return false;
}
