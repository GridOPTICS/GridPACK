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

#include <iostream>
#include <vector>

#include "base_generator_model.hpp"
#include "boost/smart_ptr/shared_ptr.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::BaseGeneratorModel::BaseGeneratorModel(void) {
  p_watch = false;
  p_hasExciter = false;
  p_hasGovernor = false;
  p_hasPss = false;
  p_hasPlantController = false;
  p_hasTorqueController = false;
  p_hasPitchController = false;
  p_hasAeroDynamicModel = false;
  p_hasDriveTrainModel = false;
  bStatus = true;
  p_generatorObservationPowerSystemBase = true;
  p_wideareafreq = 0.0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::BaseGeneratorModel::~BaseGeneratorModel(void) {}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to
 * BaseGeneratorModel
 */
void gridpack::dynamic_simulation::BaseGeneratorModel::load(
    boost::shared_ptr<gridpack::component::DataCollection> data, int idx) {}

/**
 * Update parameters in DataCollection object with current values from
 * generator
 * @param data collection object for bus that hosts generator
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::BaseGeneratorModel::updateData(
    boost::shared_ptr<gridpack::component::DataCollection> data, int idx) {}


/**
 * Initialize generator model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step
 */
void gridpack::dynamic_simulation::BaseGeneratorModel::init(double mag,
                                                            double ang,
                                                            double ts) {}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType
gridpack::dynamic_simulation::BaseGeneratorModel::INorton() {
  return gridpack::ComplexType(0.0, 0.0);
}

/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType
gridpack::dynamic_simulation::BaseGeneratorModel::NortonImpedence() {
  return gridpack::ComplexType(0.0, 0.0);
}

/**
 * Predict new state variables for time step
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseGeneratorModel::
    predictor_currentInjection(bool flag) {}

/**
 * Correct new state variables for time step
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseGeneratorModel::
    corrector_currentInjection(bool flag) {}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseGeneratorModel::predictor(double t_inc,
                                                                 bool flag) {}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::BaseGeneratorModel::corrector(double t_inc,
                                                                 bool flag) {}

void gridpack::dynamic_simulation::BaseGeneratorModel::setWideAreaFreqforPSS(
    double freq) {
  p_wideareafreq = freq;
}

bool gridpack::dynamic_simulation::BaseGeneratorModel::tripGenerator() {
  return true;
}

/**
 * return true if modify the generator parameters successfully
 * input controlTyp: 0: GFI mp adjust; 1: GFI mq adjust; 2: GFI Pset adjust; 3:
 * GFI Qset adjust; others: invalid input newParValScaletoOrg:  GFI new
 * parameter scale factor to the very initial parameter value at the begining of
 * dynamic simulation
 *
 */
bool gridpack::dynamic_simulation::BaseGeneratorModel::
    applyGeneratorParAdjustment(int controlType, double newParValScaletoOrg) {
  return false;
}

/**
 * Set voltage on each generator
 */
void gridpack::dynamic_simulation::BaseGeneratorModel::setVoltage(
    gridpack::ComplexType voltage) {}

/**
 * Set terminal voltage frequency on each generator
 */
void gridpack::dynamic_simulation::BaseGeneratorModel::setFreq(double dFreq) {}

/**
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double gridpack::dynamic_simulation::BaseGeneratorModel::getFieldVoltage() {
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
bool gridpack::dynamic_simulation::BaseGeneratorModel::serialWrite(
    char *string, const int bufsize, const char *signal) {
  string[0] = '\0';
  return false;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void gridpack::dynamic_simulation::BaseGeneratorModel::write(const char *signal,
                                                             char *string) {}

double gridpack::dynamic_simulation::BaseGeneratorModel::getAngle() {
  return 0.0;
}

void gridpack::dynamic_simulation::BaseGeneratorModel::setGovernor(
    boost::shared_ptr<BaseGovernorModel> &governor) {
  p_governor = governor;
  p_hasGovernor = true;
}

void gridpack::dynamic_simulation::BaseGeneratorModel::setExciter(
    boost::shared_ptr<BaseExciterModel> &exciter) {
  p_exciter = exciter;
  p_hasExciter = true;
}

void gridpack::dynamic_simulation::BaseGeneratorModel::setPlantController(
    boost::shared_ptr<BasePlantControllerModel> &plant) {
  p_plant = plant;
  p_hasPlantController = true;
}

void gridpack::dynamic_simulation::BaseGeneratorModel::setTorqueController(
    boost::shared_ptr<BaseMechanicalModel> &torquecontroller) {
  p_torquecontroller = torquecontroller;
  p_hasTorqueController = true;
}

void gridpack::dynamic_simulation::BaseGeneratorModel::setPitchController(
    boost::shared_ptr<BaseMechanicalModel> &pitchcontroller) {
  p_pitchcontroller = pitchcontroller;
  p_hasPitchController = true;
}

void gridpack::dynamic_simulation::BaseGeneratorModel::setDriveTrainModel(
    boost::shared_ptr<BaseMechanicalModel> &drivetrain) {
  p_drivetrainmodel = drivetrain;
  p_hasDriveTrainModel = true;
}

void gridpack::dynamic_simulation::BaseGeneratorModel::setAeroDynamicModel(
    boost::shared_ptr<BaseMechanicalModel> &aerodynamicmodel) {
  p_aerodynamicmodel = aerodynamicmodel;
  p_hasAeroDynamicModel = true;
}

void gridpack::dynamic_simulation::BaseGeneratorModel::setPss(
    boost::shared_ptr<BasePssModel> &pss) {
  p_pss = pss;
  p_hasPss = true;
}

void gridpack::dynamic_simulation::BaseGeneratorModel::AddRelay(
    boost::shared_ptr<BaseRelayModel> &p_relay) // renke add, add relay
{
  vp_relay.push_back(p_relay);
}

void gridpack::dynamic_simulation::BaseGeneratorModel::
    ClearRelay() // renke add, clear relay vector
{
  vp_relay.clear();
}

void gridpack::dynamic_simulation::BaseGeneratorModel::setWatch(bool flag) {
  p_watch = flag;
}

bool gridpack::dynamic_simulation::BaseGeneratorModel::getWatch() {
  return p_watch;
}

boost::shared_ptr<gridpack::dynamic_simulation::BaseGovernorModel>
gridpack::dynamic_simulation::BaseGeneratorModel::getGovernor() {
  return p_governor;
}

boost::shared_ptr<gridpack::dynamic_simulation::BaseExciterModel>
gridpack::dynamic_simulation::BaseGeneratorModel::getExciter() {
  return p_exciter;
}

boost::shared_ptr<gridpack::dynamic_simulation::BasePssModel>
gridpack::dynamic_simulation::BaseGeneratorModel::getPss() {
  return p_pss;
}

boost::shared_ptr<gridpack::dynamic_simulation::BasePlantControllerModel>
gridpack::dynamic_simulation::BaseGeneratorModel::getPlantController() {
  return p_plant;
}

boost::shared_ptr<gridpack::dynamic_simulation::BaseMechanicalModel>
gridpack::dynamic_simulation::BaseGeneratorModel::getTorqueController() {
  return p_torquecontroller;
}

boost::shared_ptr<gridpack::dynamic_simulation::BaseMechanicalModel>
gridpack::dynamic_simulation::BaseGeneratorModel::getPitchController() {
  return p_pitchcontroller;
}

boost::shared_ptr<gridpack::dynamic_simulation::BaseMechanicalModel>
gridpack::dynamic_simulation::BaseGeneratorModel::getDriveTrainModel() {
  return p_drivetrainmodel;
}

boost::shared_ptr<gridpack::dynamic_simulation::BaseMechanicalModel>
gridpack::dynamic_simulation::BaseGeneratorModel::getAeroDynamicModel() {
  return p_aerodynamicmodel;
}

// renke add
boost::shared_ptr<gridpack::dynamic_simulation::BaseRelayModel>
gridpack::dynamic_simulation::BaseGeneratorModel::getRelay(int iRelay) {
  return vp_relay[iRelay];
}

// renke add
void gridpack::dynamic_simulation::BaseGeneratorModel::getRelayNumber(
    int &nrelay) {
  if (vp_relay.empty()) {
    nrelay = 0;
  } else {
    nrelay = vp_relay.size();
  }
}

bool gridpack::dynamic_simulation::BaseGeneratorModel::getGenStatus() {
  return bStatus;
}

void gridpack::dynamic_simulation::BaseGeneratorModel::SetGenServiceStatus(
    bool sta) {
  bStatus = sta;
}

/**
 * return a vector containing any generator values that are being
 * watched
 * @param vals vector of watched values
 */
void gridpack::dynamic_simulation::BaseGeneratorModel::getWatchValues(
    std::vector<double> &vals) {
  vals.clear();
}

void gridpack::dynamic_simulation::BaseGeneratorModel::
    setGeneratorObPowerBaseFlag(bool generatorObservationPowerSystemBase) {
  p_generatorObservationPowerSystemBase = generatorObservationPowerSystemBase;
}

/**
 * Set internal state parameter in generator
 * @param name character string corresponding to state variable
 * @param value new value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::BaseGeneratorModel::setState(
    std::string name, double value)
{
  return false;
}

/**
 * Get internal state parameter in generator
 * @param name character string corresponding to state variable
 * @param value current value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::BaseGeneratorModel::getState(
    std::string name, double *value)
{
  return false;
}
