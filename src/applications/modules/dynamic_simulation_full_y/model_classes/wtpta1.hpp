/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wtpta1.hpp
 * @author Shuangshuang Jin 
 * @Created on:   Nov 15, 2022
 * 
 * @Last Updated: December 7, 2022
 * Shrirang Abhyankar
 * Added all the blocks
 * @brief  
 * 
 * 
 */

#ifndef _wtpta1_h_
#define _wtpta1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_mechanical_model.hpp"
#include "cblock.hpp"
#include <string>

namespace gridpack {
namespace dynamic_simulation {
class Wtpta1Model : public BaseMechanicalModel
{
  public:
    /**
     * Basic constructor
     */
    Wtpta1Model();

    /**
     * Basic destructor
     */
    virtual ~Wtpta1Model();

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
     * Set the active power  reference parameter inside the Pitch Controller model
     * @param Pref1 value of the Active power reference 
     */
    void setPref(double Pref1);
    
  /**
   * setTurbineSpeedDeviation - sets the turbine speed deviation
   * @param domega_turb : turbine speed deviation
   * From drive train model
   **/
   void setTurbineSpeedDeviation(double domega_turb);

  /**
   * setPord - sets Pord
   * @param Pord - electric Pord
   * From electrical controller model
   **/
   void setPord(double Pord);

  /**
   *  setPord0 - Sets initial power order
   *  @param Pord0 - Initial value of reference power
   *
   **/
   void setPord0(double Pord0);

  /**
   * getTheta - Get output of pitch controller
   * @output theta - output of pitch controller
   **/
   double getTheta();

  /**
   *  setOmegaref - Output of torque controller
   * @param - omega_ref : reference speed
   **/
  void setOmegaref(double omega_ref);


  private:

  // Parameters
  double Kiw, Kpw, Kic, Kpc, Kcc;
  double Tp, Thetamax, Thetamin, dThetamin, dThetamax;
  
  // Model inputs
  double Pord;
  double domega_t;
  double omega_ref;
    
  // Model output
  double Theta;
    
  // Variables
  double Pord0;
    
  // Blocks
  PIControl pitchcomp_blk;
  PIControl pitchctrl_blk;
  Filter    lag_blk;

  void computeModel(double t_inc, IntegrationStage int_flag);

};
}  // dynamic_simulation
}  // gridpack
#endif
