/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   wttqa1.cpp
 * @author Shrirang Abhyankar, Shuangshuang Jin 
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
#include "wttqa1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Wttqa1Model::Wttqa1Model(void)
{
  domega_g = 0.0;
  omega_ref = 1.0;
  Tflag = 1;
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
  if (!data->getValue(WIND_TC_TFLAG, &Tflag, idx)) Tflag = 1;  
  if (!data->getValue(WIND_TC_TP,&Tp,idx)) Tp = 0.0;
  if (!data->getValue(WIND_TC_TWREF,&Twref,idx)) Twref = 0.0; //Twref
  if (!data->getValue(WIND_TC_KPP,&Kpp,idx)) Kpp = 0.0; // Kpp
  if (!data->getValue(WIND_TC_KIP,&Kip,idx)) Kip = 0.0; // Kip
  if (!data->getValue(WIND_TC_TEMAX,&Temax,idx)) Temax = 0.0; // Temax
  if (!data->getValue(WIND_TC_TEMIN,&Temin,idx)) Temin = 0.0; //Temin
  if (!data->getValue(WIND_TC_P1,&p1,idx)) p1 = 0.0; // p1
  if (!data->getValue(WIND_TC_SPD1,&spd1,idx)) spd1 = 1.0; // spd1
  if (!data->getValue(WIND_TC_P2,&p2,idx)) p2 = 2.0; // p2
  if (!data->getValue(WIND_TC_SPD2,&spd2,idx)) spd2 = 1.0; // spd2
  if (!data->getValue(WIND_TC_P3,&p3,idx)) p3 = 0.0; // p3
  if (!data->getValue(WIND_TC_SPD3,&spd3,idx)) spd3 = 0.0; //spd3
  if (!data->getValue(WIND_TC_P4,&p4,idx)) p4 = 0.0; // p4
  if (!data->getValue(WIND_TC_SPD4,&spd4,idx)) spd4 = 0.0; // spd4
  if(!data->getValue(GENERATOR_MBASE, &MBase, idx)) MBase = 100.0;

  if(!data->getValue(WIND_TC_TRATE,&Trate,idx)) Trate = MBase;
  if(fabs(Trate) < 1e-6) Trate = MBase;

  // Set parameters in the blocks
  Pelec_filter_blk.setparams(1.0,Tp);
  wref_filter_blk.setparams(1.0,Twref);
  Tref_pi_blk.setparams(Kpp,Kip,Temin,Temax,-1000.0,1000.0);

  double x[4],y[4];
  x[0] = p1; y[0] = spd1;
  x[1] = p2; y[1] = spd2;
  x[2] = p3; y[2] = spd3;
  x[3] = p4; y[3] = spd4;

  Pomega_blk.setparams(4,x,y);
}

/**
 * Initialize  model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Wttqa1Model::init(double mag, double ang, double ts)
{
  // Initialize Pelec filter block
  Pelec_filter_blk_out = Pelec_filter_blk.init_given_u(Pelec);

  Pomega_blk_out = Pomega_blk.getoutput(Pelec_filter_blk_out);
  omega_ref = wref_filter_blk.init_given_u(Pomega_blk_out);
  // omega_ref is the same as the generator speed that we want to
  // keep the turbine at. This is the generator speed we want
  // to drive at depending on the electric power input

  // speed deviation
  domega_g = omega_ref - 1.0; // note speed deviation is negative

  Pref = Pelec;
  
  Tref_pi_blk.init_given_y(Pelec/(1+domega_g));
 
}

void gridpack::dynamic_simulation::Wttqa1Model::computeModel(double t_inc, IntegrationStage int_flag)
{
  double dtorque;
  bool   updatestate = !Vdip;
  
  Pelec_filter_blk_out = Pelec_filter_blk.getoutput(Pelec,t_inc,int_flag,true);

  Pomega_blk_out = Pomega_blk.getoutput(Pelec_filter_blk_out);
  omega_ref = wref_filter_blk.getoutput(Pomega_blk_out,t_inc,int_flag,true);

  if(Tflag == 1) {
    dtorque = (Pref0 - Pelec_filter_blk_out)/(1 + domega_g);
    Tref_pi_blk_out = Tref_pi_blk.getoutput(dtorque,t_inc,int_flag,updatestate);
  } else {
    Tref_pi_blk_out = Tref_pi_blk.getoutput(1+domega_g-omega_ref,t_inc,int_flag,updatestate);
  }

  Pref = Tref_pi_blk_out*(1 + domega_g);
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wttqa1Model::predictor(double t_inc, bool flag)
{
  computeModel(t_inc,PREDICTOR);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wttqa1Model::corrector(double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR);
}

/**
 *  Set Pref0 - Sets reference power
 *  @param Pref0 - reference power
 *  Set by plant controller model
 **/
void gridpack::dynamic_simulation::Wttqa1Model::setPref0(double Pref0_in)
{
  Pref0 = Pref0_in;
}

/**
 *  setPelec - Electrical power Pelec
 *  @param Pelec - Electrical power input
 **/
void gridpack::dynamic_simulation::Wttqa1Model::setPelec(double Pg)
{
  Pelec = Pg;
}

/**
 *  setGeneratorSpeedDeviation - Set the speed deviation
 *  @param - domega : speed deviation
 *  From drive train model
 **/
void gridpack::dynamic_simulation::Wttqa1Model::setGeneratorSpeedDeviation(double domega_in)
{
  domega_g = domega_in;
}

/**
 * setVdip - Voltage dip flag
 * @param vdip - flag to indicate voltage dip
 * From elecrical controller
 **/
void gridpack::dynamic_simulation::Wttqa1Model::setVdip(bool Vdip_in)
{
  Vdip = Vdip_in;
}

/**
 * getPref - Output of torque controller
 * @param  - Pref : reference power
 **/
double gridpack::dynamic_simulation::Wttqa1Model::getPref()
{
  return Pref;
}

/**
 *  getOmegaref - Output of torque controller
 * @param - omega_ref : reference speed
 **/
double gridpack::dynamic_simulation::Wttqa1Model::getOmegaref()
{
  return omega_ref;
}
