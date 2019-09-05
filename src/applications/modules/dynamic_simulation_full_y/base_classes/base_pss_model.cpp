/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   base_pss_model.hpp
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
#include "base_pss_model.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::BasePssModel::BasePssModel(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::BasePssModel::~BasePssModel(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * BasePssModel
 */
void gridpack::dynamic_simulation::BasePssModel::load(
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
void gridpack::dynamic_simulation::BasePssModel::init(
    double mag, double ang, double ts)
{
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BasePssModel::predictor(
    double t_inc, bool flag)
{
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BasePssModel::corrector(
    double t_inc, bool flag)
{
}

/** 
 * Get the value of the field current parameter
 * @return value of field current
 */
double gridpack::dynamic_simulation::BasePssModel::
getVstab()
{
  return 0.0;
}

void gridpack::dynamic_simulation::BasePssModel::setOmega(double omega)
{
}

double gridpack::dynamic_simulation::BasePssModel::getBusFreq(int busnum)
{
	return 0.0;
}

void gridpack::dynamic_simulation::BasePssModel::setWideAreaFreqforPSS(double freq)
{
}	
