/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_generator_model.hpp
 * @author Bruce Palmer
 * @Last modified:   May 18, 2015
 *
 * @Modified: Added wind mechanical models
 * @Last Modified: December 6, 2022
 * Shrirang Abhyankar
 *
 * @brief  Base class for generator models
 *
 *
 */

#ifndef _base_generator_model_h_
#define _base_generator_model_h_

#include "base_exciter_model.hpp"
#include "base_governor_model.hpp"
#include "base_mechanical_model.hpp"
#include "base_plant_model.hpp"
#include "base_pss_model.hpp"
#include "base_relay_model.hpp"
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"

namespace gridpack {
namespace dynamic_simulation {
class BaseGeneratorModel {
public:
  /**
   * Basic constructor
   */
  BaseGeneratorModel();

  /**
   * Basic destructor
   */
  virtual ~BaseGeneratorModel();

  /**
   * Load parameters from DataCollection object into generator model
   * @param data collection of generator parameters from input files
   * @param index of generator on bus
   */
  virtual void load(boost::shared_ptr<gridpack::component::DataCollection> data,
                    int idx);

  /**
   * Update parameters in DataCollection object with current values from
   * generator
   * @param data collection object for bus that hosts generator
   * @param index of generator on bus
   */
  virtual void updateData(boost::shared_ptr<gridpack::component::DataCollection> data,
                    int idx);

  /**
   * Initialize generator model before calculation
   * @param mag voltage magnitude
   * @param ang voltage angle
   * @param ts time step
   */
  virtual void init(double mag, double ang, double ts);

  /**
   * Return contribution to Norton current
   * @return contribution to Norton vector
   */
  virtual gridpack::ComplexType INorton();

  /**
   * Return Norton impedence
   * @return value of Norton impedence
   */
  virtual gridpack::ComplexType NortonImpedence();

  /**
   * Predict new state variables for time step
   * @param flag initial step if true
   */
  virtual void predictor_currentInjection(bool flag);

  /**
   * Correct state variables for time step
   * @param flag initial step if true
   */
  virtual void corrector_currentInjection(bool flag);

  /**
   * Predict new state variables for time step
   * @param t_inc time step increment
   * @param flag initial step if true
   */
  virtual void predictor(double t_inc, bool flag);

  /**
   * Correct state variables for time step
   * @param t_inc time step increment
   * @param flag initial step if true
   */
  virtual void corrector(double t_inc, bool flag);

  /**
   * Set voltage on each generator
   */
  virtual void setVoltage(gridpack::ComplexType voltage);

  /**
   * Set terminal voltage frequency on each generator
   */
  virtual void setFreq(double dFreq);

  /**
   * Write output from generators to a string.
   * @param string (output) string with information to be printed out
   * @param bufsize size of string buffer in bytes
   * @param signal an optional character string to signal to this
   * routine what about kind of information to write
   * @return true if bus is contributing string to output, false otherwise
   */
  virtual bool serialWrite(char *string, const int bufsize, const char *signal);

  /**
   * Get the value of the field voltage parameter
   * @return value of field voltage
   */
  virtual double getFieldVoltage();

  virtual double getAngle();
  virtual void setWideAreaFreqforPSS(double freq);

  /**
   * return true if trip generator successfully
   *
   */
  virtual bool tripGenerator();

  /**
   * return true if modify the generator parameters successfully
   * input controlTyp: 0: GFI mp adjust; 1: GFI mq adjust; 2: GFI Pset adjust;
   * 3: GFI Qset adjust; others: invalid input newParValScaletoOrg:  GFI new
   * parameter scale factor to the very initial parameter value at the begining
   * of dynamic simulation
   *
   */
  virtual bool applyGeneratorParAdjustment(int controlType,
                                           double newParValScaletoOrg);

  /**
   * Write out generator state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
  virtual void write(const char *signal, char *string);

  void setGovernor(boost::shared_ptr<BaseGovernorModel> &p_governor);

  void setExciter(boost::shared_ptr<BaseExciterModel> &p_exciter);

  void setPss(boost::shared_ptr<BasePssModel> &p_pss);

  void setPlantController(boost::shared_ptr<BasePlantControllerModel> &p_plant);

  void
  setTorqueController(boost::shared_ptr<BaseMechanicalModel> &torquecontroller);

  void
  setPitchController(boost::shared_ptr<BaseMechanicalModel> &pitchcontroller);

  void setDriveTrainModel(boost::shared_ptr<BaseMechanicalModel> &drivetrain);

  void setAeroDynamicModel(boost::shared_ptr<BaseMechanicalModel> &aerodyn);

  void
  AddRelay(boost::shared_ptr<BaseRelayModel> &p_relay); // renke add, add relay
  void ClearRelay(); // renke add, clear relay vector

  boost::shared_ptr<BaseGovernorModel> getGovernor();

  boost::shared_ptr<BaseExciterModel> getExciter();

  boost::shared_ptr<BasePssModel> getPss();

  boost::shared_ptr<BasePlantControllerModel> getPlantController();

  boost::shared_ptr<BaseMechanicalModel> getTorqueController();
  boost::shared_ptr<BaseMechanicalModel> getPitchController();
  boost::shared_ptr<BaseMechanicalModel> getDriveTrainModel();
  boost::shared_ptr<BaseMechanicalModel> getAeroDynamicModel();

  boost::shared_ptr<BaseRelayModel> getRelay(int iRelay); // renke add
  void getRelayNumber(int &nrelay);                       // renke add

  void setWatch(bool flag);

  bool getWatch();

  void setGeneratorObPowerBaseFlag(bool generatorObservationPowerSystemBase);

  /**
   * return the bolean indicating whether the gen is tripped by a relay
   */
  bool getGenStatus();

  /**
   * set the status to be false if gen is tripped by a relay
   */
  void SetGenServiceStatus(bool sta);

  /**
   * return a vector containing any generator values that are being
   * watched
   * @param vals vector of watched values
   */
  virtual void getWatchValues(std::vector<double> &vals);

  /**
   * Set internal state parameter in generator
   * @param name character string corresponding to state variable
   * @param value new value for state parameter
   * @return false if no variable corresponding to name is found
   */
  virtual bool setState(std::string name, double value);

  /**
   * Get internal state parameter in generator
   * @param name character string corresponding to state variable
   * @param value current value for state parameter
   * @return false if no variable corresponding to name is found
   */
  virtual bool getState(std::string name, double *value);


  bool p_hasExciter;
  bool p_hasGovernor;
  bool p_hasPss;
  bool p_hasPlantController;
  bool p_hasTorqueController;
  bool p_hasPitchController;
  bool p_hasDriveTrainModel;
  bool p_hasAeroDynamicModel;

  bool p_generatorObservationPowerSystemBase;
  double p_wideareafreq;

private:
  boost::shared_ptr<BaseGovernorModel> p_governor;     // governor
  boost::shared_ptr<BaseExciterModel> p_exciter;       // exciter
  boost::shared_ptr<BasePssModel> p_pss;               // stabilizer
  boost::shared_ptr<BasePlantControllerModel> p_plant; // plant controller

  boost::shared_ptr<BaseMechanicalModel>
      p_torquecontroller;                                   // Torque controller
  boost::shared_ptr<BaseMechanicalModel> p_pitchcontroller; // Pitch controller
  boost::shared_ptr<BaseMechanicalModel> p_drivetrainmodel; // Drive train model
  boost::shared_ptr<BaseMechanicalModel>
      p_aerodynamicmodel; // Aerodynamic model

  bool p_watch;
  bool bStatus;
  std::vector<boost::shared_ptr<BaseRelayModel> >
      vp_relay; // renke add, relay vector
};
} // namespace dynamic_simulation
} // namespace gridpack
#endif
