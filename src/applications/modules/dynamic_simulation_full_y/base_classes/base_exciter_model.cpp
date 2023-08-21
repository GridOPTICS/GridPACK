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
 * Set generator active and reactive power (machine MVA base)
 * 
 */
void gridpack::dynamic_simulation::BaseExciterModel::setGeneratorPower(double Pg, double Qg)
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

// Yuan added below 2020-6-23
/** 
 * Set the exciter bus number
 * @return value of exciter bus number
 */
void gridpack::dynamic_simulation::BaseExciterModel::setExtBusNum(int ExtBusNum)
{
}	

/** 
 * Set the exciter generator id
 * @return value of generator id
 */
void gridpack::dynamic_simulation::BaseExciterModel::setExtGenId(std::string ExtGenId)
{
}	

void gridpack::dynamic_simulation::BaseExciterModel::setIpcmdIqcmd(double ipcmd, double iqcmd)
{
}

void gridpack::dynamic_simulation::BaseExciterModel::setPrefQext(double pref, double qext)
{
}


double gridpack::dynamic_simulation::BaseExciterModel::getPref(){
	return 0.0;
}

double gridpack::dynamic_simulation::BaseExciterModel::getQext(){
	return 0.0;
}

double gridpack::dynamic_simulation::BaseExciterModel::getIpcmd(){
	return 0.0;
}

double gridpack::dynamic_simulation::BaseExciterModel::getIqcmd(){
	return 0.0;
}

/**
 * Set internal state parameter in exciter
 * @param name character string corresponding to state variable
 * @param value new value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::BaseExciterModel::setState( std::string name,
    double value)
{
  return false;
}

/**
 * Get internal state parameter in exciter
 * @param name character string corresponding to state variable
 * @param value current value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::BaseExciterModel::getState( std::string name,
    double *value)
{
  return false;
}

/**
 * Write output from exciter to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if governor is contributing string to output, false otherwise
 */
bool gridpack::dynamic_simulation::BaseExciterModel::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  false;
}
