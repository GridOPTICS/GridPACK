/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   base_mechanical_model.cpp
 * @author Shuangshuang Jin, Shrirang Abhyankar
 * @Created on:   Nov 15, 2022
 *
 * @Updated: Dec 5, 2022
 * Shrirang Abhyankar
 * Added Input/output methods

 * @brief Base class for wind mechanical models
 * The same base class is used for different wind mechanical components -
 * Drive train, Torque controller, Aerodynamic torque, and Pitch controller
 *
 *
 */

#include <iostream>
#include <vector>

#include "base_mechanical_model.hpp"
#include "boost/smart_ptr/shared_ptr.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::BaseMechanicalModel::BaseMechanicalModel(void) {}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::BaseMechanicalModel::~BaseMechanicalModel(void) {}

/**
 * Load parameters from DataCollection object into mechanical model
 * @param data collection of mechanical parameters from input files
 * @param index of mechanical on bus
 * TODO: might want to move this functionality to
 * BaseMechanicalModel
 */
void gridpack::dynamic_simulation::BaseMechanicalModel::load(
    boost::shared_ptr<gridpack::component::DataCollection> data, int idx) {}

/**
 * Initialize mechanical model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step
 */
void gridpack::dynamic_simulation::BaseMechanicalModel::init(double mag,
                                                             double ang,
                                                             double ts) {}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseMechanicalModel::predictor(double t_inc,
                                                                  bool flag) {}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseMechanicalModel::corrector(double t_inc,
                                                                  bool flag) {}

/** Torque controller inputs and outputs **/

/**
 *  Set Pref0 - Sets reference power
 *  @param Pref0 - reference power
 *  Set by plant controller model
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::setPref0(double Pref0) {
}

/**
 *  setPelec - Electrical power Pelec
 *  @param Pelec - Electrical power input
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::setPelec(double Pelec) {
}

/**
 *  setPmech - Mechanical power Pmech
 *  @param Pmech - Mechanical power input
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::setPmech(double Pmech) {
}


/**
 *  setGeneratorSpeedDeviation - Set the speed deviation
 *  @param - domega : speed deviation
 *  From drive train model
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::setGeneratorSpeedDeviation(
    double domega) {}

/**
 * setVdip - Voltage dip flag
 * @param vdip - flag to indicate voltage dip
 * From elecrical controller
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::setVdip(bool Vdip) {}

/**
 * getPref - Output of torque controller
 * @param  - Pref : reference power
 **/
double gridpack::dynamic_simulation::BaseMechanicalModel::getPref() {
  return 0.0;
}

/**
 *  getOmegaref - Output of torque controller
 * @param - omega_ref : reference speed
 **/
double gridpack::dynamic_simulation::BaseMechanicalModel::getOmegaref() {
  return 1.0;
}

/**
 *  setOmegaref - Output of torque controller
 * @param - omega_ref : reference speed
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::setOmegaref(double omega_ref)
{
}


/**   Pitch controller **/

/**
 * setTurbineSpeedDeviation - sets the turbine speed deviation
 * @param domega_turb : turbine speed deviation
 * From drive train model
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::
    setTurbineSpeedDeviation(double domega_turb) {}

/**
 * setPord - sets Pord
 * @param Pord - electric Pord
 * From electrical controller model
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::setPord(double Pord) {}

/**
 *  setPord0 - Sets initial power order
 *  @param Pord0 - Initial value of reference power
 *
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::setPord0(double Pord0) {
}

/**
 * getTheta - Get output of pitch controller
 * @output theta - output of pitch controller
 **/
double gridpack::dynamic_simulation::BaseMechanicalModel::getTheta() {
  return 0.0;
}

/** Aerodynamic model **/

/**
 * setTheta - sets pitch angle
 * @param Theta - pitch angle
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::setTheta(double Theta) {
}

/**
 * getTaero - returns the aero-dynamic torque
 * @output Taero - aero dynamic torque, output of aerodynamic model
 **/
double gridpack::dynamic_simulation::BaseMechanicalModel::getTaero() {
  return 0.0;
}

/** Drive train model **/

/**
 * setTmech - sets mechanical torque
 * @param Tmech - mechanical torque
 * From aerodynamic model
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::setTmech(double Tmech) {
}

/**
 * getTmech - sets mechanical torque
 * @return Tmech - mechanical torque
 **/
double gridpack::dynamic_simulation::BaseMechanicalModel::getTmech()
{
  return 0.0;
}


/**
 * setTelec - sets electrical torque
 * @param Telec - electrical torque
 * From generator?
 **/
void gridpack::dynamic_simulation::BaseMechanicalModel::setTelec(double Tmech) {
}

/**
 * getTurbineSpeedDeviation - gets the turbine speed deviation
 * @ouput turbine speed deviation
 **/
double
gridpack::dynamic_simulation::BaseMechanicalModel::getTurbineSpeedDeviation() {
  return 0.0;
}

/**
 * getGenSpeedDeviation - gets the generator speed deviation
 * @ouput generator speed deviation
 **/
double gridpack::dynamic_simulation::BaseMechanicalModel::
    getGeneratorSpeedDeviation() {
  return 0.0;
}

/**
 * getGenRotorAngleDeviation - gets the rotor angle deviation
 * @ouput rotor angle deviation
 **/
double
gridpack::dynamic_simulation::BaseMechanicalModel::getRotorAngleDeviation() {
  return 0.0;
}
