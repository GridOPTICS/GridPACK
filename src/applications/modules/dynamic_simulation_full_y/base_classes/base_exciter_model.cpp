/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   base_exciter_model.hpp
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
#include "base_exciter_model.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::BaseExciterModel::BaseExciterModel(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::BaseExciterModel::~BaseExciterModel(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * BaseExciterModel
 */
void gridpack::dynamic_simulation::BaseExciterModel::load(
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
void gridpack::dynamic_simulation::BaseExciterModel::init(
    double mag, double ang, double ts)
{
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseExciterModel::predictor(
    double t_inc, bool flag)
{
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseExciterModel::corrector(
    double t_inc, bool flag)
{
}

/**
/**
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void gridpack::dynamic_simulation::BaseExciterModel::
setFieldVoltage(double fldv)
{
}

/**
 * Set the field current parameter inside the exciter
 * @param fldc value of the field current
 */
void gridpack::dynamic_simulation::BaseExciterModel::
setFieldCurrent(double fldc)
{
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double gridpack::dynamic_simulation::BaseExciterModel::
getFieldVoltage()
{
  return 0.0;
}

/** 
 * Get the value of the field current parameter
 * @return value of field current
 */
double gridpack::dynamic_simulation::BaseExciterModel::
getFieldCurrent()
{
  return 0.0;
}

/** 
 * Set the value of the Vterminal
 * @return value of field current
 */
void gridpack::dynamic_simulation::BaseExciterModel::setVterminal(double mag)
{
}

/** 
 * Set the value of the omega
 * @return value of field current
 */
void gridpack::dynamic_simulation::BaseExciterModel::setOmega(double omega)
{
}

/** 
 * Set the value of VComp
 * @return value of Vcomp
 */
void gridpack::dynamic_simulation::BaseExciterModel::setVcomp(double Vcomp)
{
}

void gridpack::dynamic_simulation::BaseExciterModel::setVstab(double Vstab)
{
}

void gridpack::dynamic_simulation::BaseExciterModel::setWideAreaFreqforPSS(double freq)
{
}	
