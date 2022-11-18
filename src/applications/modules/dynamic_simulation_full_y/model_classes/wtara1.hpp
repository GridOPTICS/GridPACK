/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wtara1.hpp
 * @author Shuangshuang Jin 
 * @Created on:   Nov 15, 2022
 * 
 * @Last Updated: Dec 7, 2022
 * Shrirang Abhyankar
 * Added all pieces
 *
 * @brief  : Aerodynamic model
 * 
 * 
 */

#ifndef _wtara1_h_
#define _wtara1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_mechanical_model.hpp"
#include <string>

namespace gridpack {
namespace dynamic_simulation {
class Wtara1Model : public BaseMechanicalModel
{
  public:
    /**
     * Basic constructor
     */
    Wtara1Model();

    /**
     * Basic destructor
     */
    virtual ~Wtara1Model();

    /**
     * Load parameters from DataCollection object into mechanical model
     * @param data collection of mechanical parameters from input files
     * @param index of mechanical on bus
     * TODO: might want to move this functionality to BaseMechanicalModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Initialize mechanical model before calculation
     * @param mag voltage magnitude
     * @param ang voltage angle
     * @param ts time step 
     */
    void init(double mag, double ang, double ts);

    /**
     * Predict new state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    void predictor(double t_inc, bool flag);

    /**
     * Correct state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    void corrector(double t_inc, bool flag);

  /**
   * setTheta - sets pitch angle
   * @param Theta - pitch angle
   **/
   void setTheta(double Theta);

  /**
   * setTurbineSpeedDeviation - sets the turbine speed deviation
   * @param domega_turb : turbine speed deviation
   * From drive train model
   **/
   void setTurbineSpeedDeviation(double domega_turb);

  /**
   *  setPmech - Mechanical power Pmech
   *  @param Pmech - Mechanical power input
   **/
  void setPmech(double Pmech);

  /**
   * getTaero - returns the aero-dynamic torque
   * @output Taero - aero dynamic torque, output of aerodynamic model
   **/
   double getTaero();

  private:

  // Parameters
  double Ka, Theta0;
    
  // Model inputs
  double Theta; // External pitch angle
  double domega_t; // Turbine speed deviation
  double Pmech0;

  // Model outputs
  double Taero; // Aerodynamic torque
};
}  // dynamic_simulation
}  // gridpack
#endif
