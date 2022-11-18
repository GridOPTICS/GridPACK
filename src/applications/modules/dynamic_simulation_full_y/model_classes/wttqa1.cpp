/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   wttqa1.cpp
 * @author Shuangshuang Jin 
 * @Created on:    Nov 15, 2022
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <string>
#include <cstdio>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_mechanical_model.hpp"
#include "wttqa1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Wttqa1Model::Wttqa1Model(void)
{
    Pg = 0.0;
    wt = 0.0;
    Pref0 = 0.0;
    Pref = 0.0;
    
    x0dly1 = 0.0;
    x0dly2 = 0.0;
    x0pi = 0.0;
    
    x1dly1 = 0.0;
    x1dly2 = 0.0;
    x1pi = 0.0;

    fpe = 1; // predefined?
    Tflag = 1; // predefined?
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Wttqa1Model::~Wttqa1Model(void)
{
}

/**
 * Load parameters from DataCollection object into mechanical model
 * @param data collection of mechanical parameters from input files
 * @param index of mechanical on bus
 * TODO: might want to move this functionality to
 * Wttqa1Model
 */
void gridpack::dynamic_simulation::Wttqa1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
    // missing fpe?
    if (!data->getValue(WIND_TC_TP,&Ts1,idx)) Ts1 = 0.0; //Tp
    if (!data->getValue(WIND_TC_TWREF,&Ts2,idx)) Ts2 = 0.0; //Twref
    if (!data->getValue(WIND_TC_KPP,&kp1,idx)) kp1 = 0.0; // Kpp
    if (!data->getValue(WIND_TC_KIP,&ki1,idx)) ki1 = 0.0; // Kip
    if (!data->getValue(WIND_TC_TEMAX,&Max1,idx)) Max1 = 0.0; // Temax
    if (!data->getValue(WIND_TC_TEMIN,&Min1,idx)) Min1 = 0.0; //Temin
}

/**
 * Initialize  model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Wttqa1Model::init(double mag, double ang, double ts)
{
    delayblock_dly1.init(Pg, Ts1);
    delayblock_dly1.init(wt, Ts2);
    piblock_pi.init(Pref, kp1, ki1, Max1, Min1); // out3(1): Pref
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wttqa1Model::predictor(double t_inc, bool flag)
{
    // First delay
    x0dly1 = delayblock_dly1.predictor(Pg, t_inc, flag); // Pg1 is Pg?
    
    // Convert Pe to wref using fpe
    In_d = x0dly1 * fpe;
    
    // Second delay
    x0dly2 = delayblock_dly2.predictor(In_d, t_inc, flag);
    
    // Calculate the PI input determined by the top loop
    PIin1 = wt - x0dly2; // wt1 is wt?

    // Calculate the PI input determined by the bottom loop
    PIin2 = (Pref0 - x0dly1) / wt; // wt1 is wt? Pref01 is Pref0?

    if (Tflag == 0)
        PIin = PIin1;
    else
        PIin = PIin2;

    // PI controller
    // Process PI controller
    x0pi = piblock_pi.predictor(PIin, t_inc, flag);

    Pref = x0pi;

}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wttqa1Model::corrector(double t_inc, bool flag)
{
    // First delay
    x1dly1 = delayblock_dly1.predictor(Pg, t_inc, flag); // Pg1 is Pg?
    
    // Convert Pe to wref using fpe
    In_d = x1dly1 * fpe;
    
    // Second delay
    x1dly2 = delayblock_dly2.predictor(In_d, t_inc, flag);
    
    // Calculate the PI input determined by the top loop
    PIin1 = wt - x1dly2; // wt1 is wt?

    // Calculate the PI input determined by the bottom loop
    PIin2 = (Pref0 - x1dly1) / wt; // wt1 is wt? Pref01 is Pref0?

    if (Tflag == 0)
        PIin = PIin1;
    else
        PIin = PIin2;

    // PI controller
    // Process PI controller
    x1pi = piblock_pi.predictor(PIin, t_inc, flag);

    Pref = x1pi;
}

/**
 * Set the Pg parameter inside the Torque Controller model
 * @param Pg1 value of the active power
 */
void  gridpack::dynamic_simulation::Wttqa1Model::setPg(double Pg1)
{
    Pg = Pg1;
}

/**
 * Set the wt parameter inside the Torque Controller model
 * @param wt1 value of the turbine speed
 */
void  gridpack::dynamic_simulation::Wttqa1Model::setWt(double wt1)
{
    wt = wt1;
}

/**
 * Set the Pref0 parameter inside the Torque Controller model
 * @param Pref01 value of the active power reference
 */
void  gridpack::dynamic_simulation::Wttqa1Model::setPref0(double Pref01)
{
    Pref0 = Pref01;
}

/**
 * Get the value of the active power reference
 * @return value of active power reference
 */
double  gridpack::dynamic_simulation::Wttqa1Model::getPref()
{
  return Pref;
}
