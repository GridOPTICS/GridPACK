/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   base_load_model.hpp
 * @author Bruce Palmer
 * @Last modified:   September 23, 2016
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_load_model.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::BaseLoadModel::BaseLoadModel(void)
{
  p_watch = false;
  dyn_p = 0.0;
  dyn_q = 0.0;
  dyn_load_id = " ";
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::BaseLoadModel::~BaseLoadModel(void)
{
}

/**
 * Load parameters from DataCollection object into load model
 * @param data collection of load parameters from input files
 * @param index of load on bus
 */
void gridpack::dynamic_simulation::BaseLoadModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx, double loadP, double loadQ, int ibCMPL)
{
}

/**
 * Initialize load model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step
 */
void gridpack::dynamic_simulation::BaseLoadModel::init(
    double mag, double ang, double ts) 
{
}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::BaseLoadModel::INorton()
{
  return gridpack::ComplexType(0.0,0.0);
}

/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::BaseLoadModel::NortonImpedence()
{
  return gridpack::ComplexType(0.0,0.0);
}

/**
 * Predict new state variables for time step
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseLoadModel::predictor_currentInjection(bool flag)
{
}

/**
 * Correct new state variables for time step
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseLoadModel::corrector_currentInjection(bool flag)
{
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseLoadModel::predictor(
    double t_inc, bool flag)
{
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseLoadModel::corrector(
    double t_inc, bool flag)
{
}

/**
* post process for time step
* @param t_inc time step increment
* @param flag initial step if true
*/
void gridpack::dynamic_simulation::BaseLoadModel::dynamicload_post_process(
    double t_inc, bool flag)
{
}

/**
 * Set voltage on each load
 */
void gridpack::dynamic_simulation::BaseLoadModel::setVoltage(
    gridpack::ComplexType voltage)
{
  
}

/**
 * Set terminal voltage frequency on each load
 */
void gridpack::dynamic_simulation::BaseLoadModel::setFreq(double dFreq)
{
  
}

/**
 * get intialized reactive power of the dynamic load model
 */
double gridpack::dynamic_simulation::BaseLoadModel::getInitReactivePower(void)
{
	return 0.0;
}

/**
 * Write output from loads to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::dynamic_simulation::BaseLoadModel::
serialWrite(char *string, const int bufsize, const char *signal)
{
  return false;
}

/**
 * Write out load state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void gridpack::dynamic_simulation::BaseLoadModel::write(
    const char* signal, char *string)
{
}

/**
 * Set watch parameter
 * @param flag set watch to true or false
 */
void gridpack::dynamic_simulation::BaseLoadModel::setWatch(bool flag)
{
  p_watch = flag;
}

/**
 * Get current value of watch parameter
 * @return current value of watch parameter
 */
bool gridpack::dynamic_simulation::BaseLoadModel::getWatch()
{
  return p_watch;
}

void gridpack::dynamic_simulation::BaseLoadModel::setDynLoadP(double pl)
{
	dyn_p = pl;
}
void gridpack::dynamic_simulation::BaseLoadModel::setDynLoadQ(double ql)
{
	dyn_q = ql;
}
void gridpack::dynamic_simulation::BaseLoadModel::setDynLoadID(std::string load_id)
{
	dyn_load_id = load_id;
}

double gridpack::dynamic_simulation::BaseLoadModel::getDynLoadP()
{
	return dyn_p;
}
double gridpack::dynamic_simulation::BaseLoadModel::getDynLoadQ()
{
	return dyn_q;
}
std::string gridpack::dynamic_simulation::BaseLoadModel::getDynLoadID()
{
	return dyn_load_id;
}
