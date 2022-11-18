/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wttqa1.hpp
 * @author Shrirang Abhyankar, Shuangshuang Jin
 * @Created on:   Nov 15, 2022
 * 
 * @brief  
 * 
 * 
 */

#ifndef _wttqa1_h_
#define _wttqa1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_mechanical_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"
#include <string>

namespace gridpack {
namespace dynamic_simulation {
class Wttqa1Model : public BaseMechanicalModel
{
  public:
    /**
     * Basic constructor
     */
    Wttqa1Model();

    /**
     * Basic destructor
     */
    virtual ~Wttqa1Model();

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
   *  Set Pref0 - Sets reference power
   *  @param Pref0 - reference power
   *  Set by plant controller model
   **/
   void setPref0(double Pref0);

  /**
   *  setPelec - Electrical power Pelec
   *  @param Pelec - Electrical power input
   **/
   void setPelec(double Pelec);

  /**
   *  setGeneratorSpeedDeviation - Set the speed deviation
   *  @param - domega : speed deviation
   *  From drive train model
   **/
   void setGeneratorSpeedDeviation(double domega);

  /**
   * setVdip - Voltage dip flag
   * @param vdip - flag to indicate voltage dip
   * From elecrical controller
   **/
   void setVdip(bool Vdip);
  
  /**
   * getPref - Output of torque controller
   * @param  - Pref : reference power
   **/
   double getPref();

  /**
   *  getOmegaref - Output of torque controller
   * @param - omega_ref : reference speed
   **/
   double getOmegaref();
  
  private:

  // Parameters
  int    Tflag;
  double Kpp, Kip, Tp, Twref, Temax, Temin;
  double p1, spd1, p2, spd2, p3, spd3, p4, spd4;
  double Trate;
  
  // Inputs:
  double Pelec; // Generator real power ouput
  double Pref0; // Plant controller reference
  double domega_g; // speed deviation
  bool   Vdip;   // Voltage dip flag
  
  // Outputs
  double Pref; // Reference power output
  double omega_ref; // Reference speed (pu)

  // Internal variables
  double MBase; // Machine base

  // Control blocks
  Filter Pelec_filter_blk; // Pelec filter block
  double Pelec_filter_blk_out; // output of Pelec filter block
  
  Filter wref_filter_blk; // wref filter block

  PIControl Tref_pi_blk; // PI controller for Tref
  double    Tref_pi_blk_out; // Output of Tref PI control block

  PiecewiseSlope Pomega_blk; // Power speed piecewise linear block
  double         Pomega_blk_out; // Output of Pomega_blk

  void computeModel(double t_inc,IntegrationStage flag);
};
}  // dynamic_simulation
}  // gridpack
#endif

