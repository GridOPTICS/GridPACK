/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_mechanical_model.hpp
 * @author Shuangshuang Jin, Shrirang Abhyankar
 * @Created on:    Nov 15, 2022
 *
 * @brief
 *
 *
 */

#ifndef _base_mechanical_model_h_
#define _base_mechanical_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"

namespace gridpack {
namespace dynamic_simulation {
class BaseMechanicalModel {
public:
  /**
   * Basic constructor
   */
  BaseMechanicalModel();

  /**
   * Basic destructor
   */
  virtual ~BaseMechanicalModel();

  /**
   * Load parameters from DataCollection object into mechanical model
   * @param data collection of mechanical parameters from input files
   * @param index of mechanical on bus
   * TODO: might want to move this functionality to BaseMechanicalModel
   */
  virtual void load(boost::shared_ptr<gridpack::component::DataCollection> data,
                    int idx);

  /**
   * Initialize mechanical model before calculation
   * @param mag voltage magnitude
   * @param ang voltage angle
   * @param ts time step
   */
  virtual void init(double mag, double ang, double ts);

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

  /** Torque controller inputs and outputs **/

  /**
   *  Set Pref0 - Sets reference power
   *  @param Pref0 - reference power
   *  Set by plant controller model
   **/
  virtual void setPref0(double Pref0);

  /**
   *  setPelec - Electrical power Pelec
   *  @param Pelec - Electrical power input
   **/
  virtual void setPelec(double Pelec);

  /**
   *  setPmech - Mechanical power Pmech
   *  @param Pmech - Mechanical power input
   **/
  virtual void setPmech(double Pmech);


  /**
   *  setGeneratorSpeedDeviation - Set the speed deviation
   *  @param - domega : speed deviation
   *  From drive train model
   **/
  virtual void setGeneratorSpeedDeviation(double domega);

  /**
   * setVdip - Voltage dip flag
   * @param vdip - flag to indicate voltage dip
   * From elecrical controller
   **/
  virtual void setVdip(bool Vdip);

  /**
   * getPref - Output of torque controller
   * @param  - Pref : reference power
   **/
  virtual double getPref();

  /**
   *  getOmegaref - Output of torque controller
   * @return - omega_ref : reference speed
   **/
  virtual double getOmegaref();

  /**
   *  setOmegaref - Output of torque controller
   * @param - omega_ref : reference speed
   **/
  virtual void setOmegaref(double omega_ref);


  /**   Pitch controller **/

  /**
   * setTurbineSpeedDeviation - sets the turbine speed deviation
   * @param domega_turb : turbine speed deviation
   * From drive train model
   **/
  virtual void setTurbineSpeedDeviation(double domega_turb);

  /**
   * setPord - sets Pord
   * @param Pord - electric Pord
   * From electrical controller model
   **/
  virtual void setPord(double Pord);

  /**
   *  setPord0 - Sets initial power order
   *  @param Pord0 - Initial value of reference power
   *
   **/
  virtual void setPord0(double Pord0);

  /**
   * getTheta - Get output of pitch controller
   * @output theta - output of pitch controller
   **/
  virtual double getTheta();

  /** Aerodynamic model **/

  /**
   * setTheta - sets pitch angle
   * @param Theta - pitch angle
   **/
  virtual void setTheta(double Theta);

  /**
   * getTaero - returns the aero-dynamic torque
   * @output Taero - aero dynamic torque, output of aerodynamic model
   **/
  virtual double getTaero();

  /** Drive train model **/

  /**
   * setTmech - sets mechanical torque
   * @param Tmech - mechanical torque
   * From aerodynamic model
   **/
  virtual void setTmech(double Tmech);

  /**
   * getTmech - sets mechanical torque
   * @return Tmech - mechanical torque
   **/
  virtual double getTmech();

  /**
   * setTelec - sets electrical torque
   * @param Telec - electrical torque
   * From generator?
   **/
  virtual void setTelec(double Tmech);

  /**
   * getTurbineSpeedDeviation - gets the turbine speed deviation
   * @ouput turbine speed deviation
   **/
  virtual double getTurbineSpeedDeviation();

  /**
   * getGenSpeedDeviation - gets the generator speed deviation
   * @ouput generator speed deviation
   **/
  virtual double getGeneratorSpeedDeviation();

  /**
   * getGenRotorAngleDeviation - gets the rotor angle deviation
   * @ouput rotor angle deviation
   **/
  virtual double getRotorAngleDeviation();

private:
};
} // namespace dynamic_simulation
} // namespace gridpack
#endif
