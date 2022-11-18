/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wttqa1.hpp
 * @author Shuangshuang Jin 
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
#include "PIBlockwithLimit.hpp"
#include "DelayBlock.hpp"
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
     * Set the Pg parameter inside the Torque Controller model
     * @param Pg1 value of the active power
     */
    void setPg(double Pg1);
    
    /**
     * Set the wt parameter inside the Torque Controller model
     * @param wt1 value of the turbine speed
     */
    void setWt(double wt1);
    
    /**
     * Set the Pref0 parameter inside the Torque Controller model
     * @param Pref01 value of the active power reference
     */
    void setPref0(double Pref01);
    
    /**
     * Get the value of the active power reference
     * @return value of  active power reference
     */
    double getPref();

  private:

    // Mechanical Wttqa1 Parameters read from dyr:
    double Ts1, Ts2, kp1, ki1, Max1, Min1;
    
    // Inputs:
    double Pg, wt, Pref0;

    // Outputs: Active power reference
    double Pref;
    
    // WTTQA1 state variables
    double x0dly1, x0dly2, x0pi;
    double x1dly1, x1dly2, x1pi;
    double dx0dly1, dx0dly2, dx0pi;
    double dx1dly1, dx1dly2, dx1pi;

    double In_d, PIin1, PIin2, PIin, fpe, Tflag;
    
    gridpack::dynamic_simulation::DelayBlock delayblock_dly1;
    gridpack::dynamic_simulation::DelayBlock delayblock_dly2;
    gridpack::dynamic_simulation::PIBlockwithLimit piblock_pi;

};
}  // dynamic_simulation
}  // gridpack
#endif

