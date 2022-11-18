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
 * @Last Updated: Dec 7, 2022
 * Shrirang Abhyankar
 * Added all pieces
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
#include "wtpta1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Wtpta1Model::Wtpta1Model(void)
{
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
    if (!data->getValue(WIND_PC_KIW, &Kiw, idx)) Kiw = 0.0; 
    if (!data->getValue(WIND_PC_KPW, &Kpw, idx)) Kpw = 0.0; 
    if (!data->getValue(WIND_PC_KIC, &Kic, idx)) Kic = 0.0; 
    if (!data->getValue(WIND_PC_KPC ,&Kpc,idx)) Kpc = 0.0; 
    if (!data->getValue(WIND_PC_THETAMAX,&Thetamax,idx))  Thetamax = 0.0; 
    if (!data->getValue(WIND_PC_THETAMIN,&Thetamin,idx)) Thetamin = 0.0;
    if (!data->getValue(WIND_PC_RTHETAMAX,&dThetamax,idx))  dThetamax = 0.0; 
    if (!data->getValue(WIND_PC_RTHETAMIN,&dThetamin,idx)) dThetamin = 0.0;
    if (!data->getValue(WIND_PC_TP,&Tp,idx)) Tp = 0.0;

    // Set parameters for blocks
    pitchcomp_blk.setparams(Kpc,Kic,Thetamin,Thetamax,-1000.0,1000.0);
    pitchctrl_blk.setparams(Kpw,Kiw,Thetamin,Thetamax,-1000.0,1000.0);
    lag_blk.setparams(1.0,Tp,Thetamin,Thetamax,dThetamin,dThetamax,-1000.0,1000.0);
}

/**
 * Initialize  model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Wtpta1Model::init(double mag, double ang, double ts)
{
  double y1,y2;

  y1 = pitchcomp_blk.init_given_u(Pord-Pord0);
  y2 = pitchctrl_blk.init_given_u((1+domega_t)-omega_ref + Kcc*(Pord - Pord0));

  Theta = lag_blk.init_given_u(y1+y2);
}

void gridpack::dynamic_simulation::Wtpta1Model::computeModel(double t_inc,IntegrationStage int_flag)
{
  double y1, y2;

  y1 = pitchcomp_blk.getoutput(Pord-Pord0,t_inc,int_flag,true);
  y2 = pitchctrl_blk.getoutput((1+domega_t)-omega_ref + Kcc*(Pord - Pord0),t_inc,int_flag,true);

  Theta = lag_blk.getoutput(y1+y2,t_inc,int_flag,true);
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtpta1Model::predictor(double t_inc, bool flag)
{
  computeModel(t_inc,PREDICTOR);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtpta1Model::corrector(double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR);
}

/**
 * setTurbineSpeedDeviation - sets the turbine speed deviation
 * @param domega_turb : turbine speed deviation
 * From drive train model
 **/
void gridpack::dynamic_simulation::Wtpta1Model::setTurbineSpeedDeviation(double domega_turb)
{
  domega_t = domega_turb;
}

/**
 * setPord - sets Pord
 * @param Pord - electric Pord
 * From electrical controller model
 **/
void gridpack::dynamic_simulation::Wtpta1Model::setPord(double Pord_in)
{
  Pord = Pord_in;
}

/**
 *  setPord0 - Sets initial power order
 *  @param Pord0 - Initial value of reference power
 *
 **/
void gridpack::dynamic_simulation::Wtpta1Model::setPord0(double Pord0_in)
{
  Pord0 = Pord0_in;
}

/**
 *  setOmegaref - Output of torque controller
 * @param - omega_ref : reference speed
 **/
void gridpack::dynamic_simulation::Wtpta1Model::setOmegaref(double omega_ref_in)
{
  omega_ref = omega_ref_in;
}


/**
 * getTheta - Get output of pitch controller
 * @output Theta - pitch angle (degrees) output of pitch controller
 **/
double gridpack::dynamic_simulation::Wtpta1Model::getTheta()
{
  return Theta;
}
