/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wtpta1_cblock.hpp
 * @author Shuangshuang Jin 
 * @Created on:   DEC 05, 2022
 * 
 * @brief  
 * 
 * 
 */

#ifndef _wtpta1_h_
#define _wtpta1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_mechanical_model.hpp"
//#include "PIBlockwithLimit.hpp"
//#include "DelayBlockwithLimit.hpp"
#include "cblock.hpp"
#include "dblock.hpp"
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
     * Set the turbine speed parameter inside the Pitch Controller model
     * @param wt1 values of the speed
     */
    void setSpeed(double wt1);
    
    /**
     * Set the turbine speed reference parameter inside the Pitch Controller model
     * @param wref1 values of the speed reference
     */
    void setSpeedReference(double wref1);
    
    /**
     * Set the active power order parameter inside the Pitch Controller model
     * @param Pord1 value of the Active power order
     */
    void setPord(double Pord1);
    
    /**
     * Set the active power  reference parameter inside the Pitch Controller model
     * @param Pref1 value of the Active power reference 
     */
    void setPref(double Pref1);
    
    /**
     * Get the value of the pitch angle
     * @return value of  pitch angle
     */
    double getPitchAngle();

  private:

    // Mechanical Wtpta1 Parameters read from dyr:
    double Kcc, Kiw, Kpw, Kic, Kpc, Max1, Min1, RMax1, RMin1, Ts;
    
    // Input: external turbine speed, turbine speed reference, active power order
    double wt, wref, Pord;
    double Pref0; // external active power reference

    // Outputs: Pitch angle
    double Theta;
    
    // WTPTA1 state variables
    /*double x0piup, x0pidown, x0dly;
    double x1piup, x1pidown, x1dly;
    double dx0piup, dx0pidown, dx0dly;
    double dx1piup, dx1pidown, dx1dly;*/
    
    double In_d;
    
    /*gridpack::dynamic_simulation::PIBlockwithLimit piblock_up;
    gridpack::dynamic_simulation::PIBlockwithLimit piblock_down;
    gridpack::dynamic_simulation::DelayBlockwithLimit delayblock_dly;*/
    PIControl piblock_up, piblock_down;
    Filter delayblock_dly;

};
}  // dynamic_simulation
}  // gridpack
#endif
