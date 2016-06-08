/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   base_generator_model.hpp
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
#include "gridpack/include/gridpack.hpp"
#include "base_generator_model.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::DSFBaseGeneratorModel::DSFBaseGeneratorModel(void)
{
  p_watch = false;
  p_hasExciter = false;
  p_hasGovernor = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::DSFBaseGeneratorModel::~DSFBaseGeneratorModel(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to
 * DSFBaseGeneratorModel
 */
void gridpack::dynamic_simulation::DSFBaseGeneratorModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
}

/**
 * Initialize generator model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step
 */
void gridpack::dynamic_simulation::DSFBaseGeneratorModel::init(
    double mag, double ang, double ts) 
{
}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::DSFBaseGeneratorModel::INorton()
{
  return gridpack::ComplexType(0.0,0.0);
}

/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::DSFBaseGeneratorModel::NortonImpedence()
{
  return gridpack::ComplexType(0.0,0.0);
}

/**
 * Predict new state variables for time step
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFBaseGeneratorModel::predictor_currentInjection(bool flag)
{
}

/**
 * Correct new state variables for time step
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFBaseGeneratorModel::corrector_currentInjection(bool flag)
{
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFBaseGeneratorModel::predictor(
    double t_inc, bool flag)
{
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFBaseGeneratorModel::corrector(
    double t_inc, bool flag)
{
}

/**
 * Set voltage on each generator
 */
void gridpack::dynamic_simulation::DSFBaseGeneratorModel::setVoltage(
    gridpack::ComplexType voltage)
{
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double gridpack::dynamic_simulation::DSFBaseGeneratorModel::
getFieldVoltage()
{
  return 0.0;
}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::dynamic_simulation::DSFBaseGeneratorModel::
serialWrite(char *string, const int bufsize, const char *signal)
{
  return false;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void gridpack::dynamic_simulation::DSFBaseGeneratorModel::write(
    const char* signal, char *string)
{
}

double gridpack::dynamic_simulation::DSFBaseGeneratorModel::
getAngle()
{
  return 0.0;
}

void
gridpack::dynamic_simulation::DSFBaseGeneratorModel::setGovernor(boost::shared_ptr<DSFBaseGovernorModel>
    &governor)
{
  p_governor = governor;
}

void
gridpack::dynamic_simulation::DSFBaseGeneratorModel::setExciter(boost::shared_ptr<DSFBaseExciterModel>
    &exciter)
{
  p_exciter = exciter;
}

void gridpack::dynamic_simulation::DSFBaseGeneratorModel::setWatch(bool flag)
{
  p_watch = flag;
}

bool gridpack::dynamic_simulation::DSFBaseGeneratorModel::getWatch()
{
  return p_watch;
}

boost::shared_ptr<gridpack::dynamic_simulation::DSFBaseGovernorModel>
gridpack::dynamic_simulation::DSFBaseGeneratorModel::getGovernor()
{
  return p_governor;
}

boost::shared_ptr<gridpack::dynamic_simulation::DSFBaseExciterModel>
gridpack::dynamic_simulation::DSFBaseGeneratorModel::getExciter()
{
  return p_exciter;
}


