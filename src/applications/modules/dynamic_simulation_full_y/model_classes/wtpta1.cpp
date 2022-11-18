/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   wtpta1.cpp
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
#include "wtpta1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Wtpta1Model::Wtpta1Model(void)
{
    wt = 0.0;
    wref = 0.0;
    Pord = 0.0;
    Pref0 = 0.0;
    In_d = 0.0;
    
    x0piup = 0.0;
    x0pidown = 0.0;
    x0dly = 0.0;
    
    x1piup = 0.0;
    x1pidown = 0.0;
    x1dly = 0.0;
    
    Kcc = 0;
    
    Theta = 0.0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Wtpta1Model::~Wtpta1Model(void)
{
}

/**
 * Load parameters from DataCollection object into mechanical model
 * @param data collection of mechanical parameters from input files
 * @param index of mechanical on bus
 * TODO: might want to move this functionality to
 * Wtpta1Model
 */
void gridpack::dynamic_simulation::Wtpta1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
    if (!data->getValue(WIND_PC_KCC, &Kcc, idx)) Kcc = 0.0;
    if (!data->getValue(WIND_PC_KIW, &Kiw, idx)) Kiw = 0.0; //kiup
    if (!data->getValue(WIND_PC_KPW, &Kpw, idx)) Kpw = 0.0; //kpup
    if (!data->getValue(WIND_PC_KIC, &Kic, idx)) Kic = 0.0; //kido
    if (!data->getValue(WIND_PC_KPC ,&Kpc,idx)) Kpc = 0.0; //kpdo
    if (!data->getValue(WIND_PC_THETAMAX,&Max1,idx)) Max1 = 0.0; //pimax
    if (!data->getValue(WIND_PC_THETAMIN,&Min1,idx)) Min1 = 0.0; //pimin
    if (!data->getValue(WIND_PC_TP,&Ts,idx)) Ts = 0.0; //Tp, Tpi
}

/**
 * Initialize  model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Wtpta1Model::init(double mag, double ang, double ts)
{
    //piblock_up.init(wt - wref, kpup, kiup, Max1, Min1); // In1 = wt - wref ?
    piblock_up.init(wt - wref, Kpw, Kiw, Max1, Min1); // In1 = wt - wref ?
    //piblock_down.init(Pord - Pref0, kpdo, kido, Max1, Min1); //In2 = Pord - Pref0 ?
    piblock_down.init(Pord - Pref0, Kpc, Kic, Max1, Min1); //In2 = Pord - Pref0 ?
    delayblock_dly.init(Theta, Ts, Max1, Min1);
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtpta1Model::predictor(double t_inc, bool flag)
{
    x0piup = piblock_up.predictor(wt - wref, t_inc, flag);
    x0pidown = piblock_down.predictor(Pord - Pref0, t_inc, flag);

    In_d = x0piup + x0pidown;

    x0dly = delayblock_dly.predictor(In_d, t_inc, flag);

    Theta = x0dly;
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtpta1Model::corrector(double t_inc, bool flag)
{
    x1piup = piblock_up.predictor(wt - wref, t_inc, flag);
    x1pidown = piblock_down.predictor(Pord - Pref0, t_inc, flag);

    In_d = x1piup + x1pidown;

    x1dly = delayblock_dly.predictor(In_d, t_inc, flag);

    Theta = x1dly;
}

/**
 * Set the turbine speed parameter inside the Pitch Controller model
 * * @param wt1 values of the speed
 */
void  gridpack::dynamic_simulation::Wtpta1Model::setSpeed(double wt1)
{
    wt = wt1;
}

/**
 * Set the turbine speed reference parameter inside the Pitch Controller model
 * @param wref1 values of the speed reference
 */
void  gridpack::dynamic_simulation::Wtpta1Model::setSpeedReference(double wref1)
{
    wref = wref1;
}

/**
 * Set the active power order parameter inside the Pitch Controller model
 * @param Pord1 value of the Active power order
 */
void  gridpack::dynamic_simulation::Wtpta1Model::setPord(double Pord1)
{
    Pord = Pord1;
}

/**
 * Set the active power  reference parameter inside the Pitch Controller model
 * @param Pref1 value of the Active power reference 
 */
void  gridpack::dynamic_simulation::Wtpta1Model::setPref(double Pref1)
{
    Pref0 = Pref1;
}

/**
 * Get the value of the pitch angle
 * @return value of pitch angle
 */
double  gridpack::dynamic_simulation::Wtpta1Model::getPitchAngle()
{
  return Theta;
}
