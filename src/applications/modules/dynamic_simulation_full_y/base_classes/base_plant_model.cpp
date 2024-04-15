/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   base_plant_model.hpp
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
#include "base_plant_model.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::BasePlantControllerModel::BasePlantControllerModel(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::BasePlantControllerModel::~BasePlantControllerModel(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * BasePlantControllersModel
 */
void gridpack::dynamic_simulation::BasePlantControllerModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::BasePlantControllerModel::init(
    double mag, double ang, double ts)
{
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BasePlantControllerModel::predictor(
    double t_inc, bool flag)
{
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BasePlantControllerModel::corrector(
    double t_inc, bool flag)
{
}


void gridpack::dynamic_simulation::BasePlantControllerModel::setGenPQV(double P, double Q, double Vt)
{
}

void gridpack::dynamic_simulation::BasePlantControllerModel::setBusFreq(double freq)
{
}

void gridpack::dynamic_simulation::BasePlantControllerModel::setPrefQext(double Pref, double Qext)
{
}

void gridpack::dynamic_simulation::BasePlantControllerModel::setExtBusNum(int ExtBusNum)
{
}

void gridpack::dynamic_simulation::BasePlantControllerModel::setExtGenId(std::string ExtGenId)
{
}

double gridpack::dynamic_simulation::BasePlantControllerModel::getPref( )
{
	return 0.0;
}

double gridpack::dynamic_simulation::BasePlantControllerModel::getQext( )
{
	return 0.0;
}

/**
 * Set internal state parameter in plant controller
 * @param name character string corresponding to state variable
 * @param value new value for state parameter
 * @return false if no variable corresponding to name is
 found
 */
bool gridpack::dynamic_simulation::BasePlantControllerModel::setState( std::string name,
    double value)
{
  return false;
}

/**
 * Get internal state parameter in plant controller
 * @param name character string corresponding to state variable
 * @param value current value for state parameter
 * @return false if no variable corresponding to name is
 found
 */
bool gridpack::dynamic_simulation::BasePlantControllerModel::getState( std::string name,
    double *value)
{
  return false;
}

/**
 * Write output from plant model to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if governor is contributing string to output, false otherwise
 */
bool gridpack::dynamic_simulation::BasePlantControllerModel::serialWrite(
    char *string, const int bufsize, const char *signal)
{
  false;
}
