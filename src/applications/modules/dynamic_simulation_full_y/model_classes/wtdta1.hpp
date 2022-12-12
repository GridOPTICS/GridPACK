/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wtdta1.hpp
 * @author Shuangshuang Jin 
 * @Created on:   Nov 15, 2022
 * 
 * @Updated
 * Shrirang Abhyankar
 * Added in all different pieces
 * @brief  Wind drive train model
 * 
 * 
 */

#ifndef _wtdta1_h_
#define _wtdta1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_mechanical_model.hpp"
#include "cblock.hpp"
#include <string>

namespace gridpack {
namespace dynamic_simulation {
class Wtdta1Model : public BaseMechanicalModel
{
  public:
    /**
     * Basic constructor
     */
    Wtdta1Model();

    /**
     * Basic destructor
     */
    virtual ~Wtdta1Model();

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
   * setTmech - sets mechanical torque
   * @param Tmech - mechanical torque
   * From aerodynamic model
   **/
   void setTmech(double Tmech);
  
  /**
   * setTelec - sets electrical torque
   * @param Telec - electrical torque
   * From generator?
   **/
   void setTelec(double Tmech);

  /**
   * getTurbineSpeedDeviation - gets the turbine speed deviation
   * @ouput turbine speed deviation
   **/
   double getTurbineSpeedDeviation();

  /**
   * getGenSpeedDeviation - gets the generator speed deviation
   * @ouput generator speed deviation
   **/
   double getGeneratorSpeedDeviation();

  /**
   * getGenRotorAngleDeviation - gets the rotor angle deviation
   * @ouput rotor angle deviation
   **/
   double getRotorAngleDeviation();

  
  private:

  // Parameters
  double H, D, Hfrac;
  double Freq1, Dshaft;
  
  // Model inputs
  double Tm; // Mechanical torque
  double Te; // Electrical torque
  double s0; // Initial slip

  // Model outputs
  double domega_g; // speed deviation
  double dtheta_g; // rotor angle speed deviation
  double domega_t; // Turbine speed deviation

  // Blocks
  Integrator domegag_blk;
  Integrator domegat_blk;
  Integrator dthetag_blk;
  Integrator Tshaft_blk;

  // Internal variables
  double Hg; // generator inertia
  double Ht; // turbine inertia
  double Kshaft;
  bool   double_mass; // True -> also includes turbine piece

  /**
     computeModel - Used both by predictor and corrector
  **/
  void computeModel(double t_inc,IntegrationStage int_flag);
};
}  // dynamic_simulation
}  // gridpack
#endif
