/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   ieeet1.cpp
 * @author Shrirang Abhyankar
 * @Updated:   December 25, 2022
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <cstdio>

#include <cstring>
#include <string>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_exciter_model.hpp"
#include "ieeet1.hpp"

#define TS_THRESHOLD 4

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Ieeet1Model::Ieeet1Model(void)
{
  Vs = 0.0;
  sat_A = 0.0;
  sat_B = 0.0;
  has_Sat = true;
  zero_TA = false;
  zero_TR = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Ieeet1Model::~Ieeet1Model(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 */
void gridpack::dynamic_simulation::Ieeet1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(EXCITER_TR, &TR, idx)) TR = 0.0; // TR
  if (!data->getValue(EXCITER_KA, &KA, idx)) KA = 0.0; // KA 
  if (!data->getValue(EXCITER_TA, &TA, idx)) TA = 0.0; // TA
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 0.0; // Vrmax
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = 0.0; // Vrmin
  if (!data->getValue(EXCITER_KE, &KE, idx)) KE = 0.0; // KE
  if (!data->getValue(EXCITER_TE, &TE, idx)) TE = 0.0; // TE
  if (!data->getValue(EXCITER_KF, &KF, idx)) KF = 0.0; // KF
  if (!data->getValue(EXCITER_TF1, &TF1, idx)) TF1 = 0.0; // TF
  if (!data->getValue(EXCITER_E1, &E1, idx)) E1 = 0.0; // E1
  if (!data->getValue(EXCITER_SE1, &SE1, idx)) SE1 = 0.0; // SE1
  if (!data->getValue(EXCITER_E2, &E2, idx)) E2 = 0.0; // E2
  if (!data->getValue(EXCITER_SE2, &SE2, idx)) SE2 = 0.0; // SE2

  if(fabs(SE1*SE2) < 1e-6) has_Sat = false;
  if(fabs(TA) < 1e-6) zero_TA = true;
  if(fabs(TR) < 1e-6) zero_TR = true;

  if(has_Sat) {
    /* Calculate saturation function constants */
    double r;
    r = std::sqrt(SE1*E1/(SE2*E2));
    sat_A = (E1 - r*E2)/(1 - r);
    sat_B = SE1*E1/((E1 - sat_A)*(E1 - sat_A));
  }

  // Set up blocks
  if(!zero_TR) {
    Vmeas_blk.setparams(1.0,TR);
  }

  if(!zero_TA) {
    Regulator_blk.setparams(KA,TA,Vrmin,Vrmax,-1000.0,1000.0);
  } else {
    Regulator_gain_blk.setparams(KA,Vrmin,Vrmax);
  }

  double a[2],b[2];
  a[0] = TF1; a[1] = 1.0;
  b[0] = KF;  b[1] = 0.0;
  Feedback_blk.setcoeffs(a,b);

  Output_blk.setparams(TE);
}

/**
 *  * Saturation function
 *  * @ param x
 *  * return - saturation value given by quadratic function 
 *       S = B(x - A)^2/x if x > A
 *  else S = 0
 *    */
double gridpack::dynamic_simulation::Ieeet1Model::Sat(double x)
{
  double S;

  if(x < sat_A) S = 0.0;
  else S = sat_B*(x - sat_A)*(x - sat_A)/x;

  return S;
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Ieeet1Model::init(double Vm, double Va, double ts)
{
  VF = Feedback_blk.init_given_u(Efd);

  double output_blk_in;
  output_blk_in = Output_blk.init_given_y(Efd);

  double sat_signal = KE*Efd;
  if(has_Sat) {
    sat_signal += Efd*Sat(Efd);
  }

  VR = output_blk_in + sat_signal;

  double Regulator_blk_in;
  if(zero_TA) {
    Regulator_blk_in = VR/KA;
  } else {
    Regulator_blk_in = Regulator_blk.init_given_y(VR);
  }

  double Verr;
  Verr = Regulator_blk_in + VF;
      
  if(!zero_TR) {
    Vmeas = Vmeas_blk.init_given_u(Ec);
  } else Vmeas = Ec;

  Vref = Verr + Vmeas - Vs;
}

/**
 * computeModel - Updates the model states and gets the output
 */
void gridpack::dynamic_simulation::Ieeet1Model::computeModel(double t_inc,IntegrationStage int_flag)
{
  if(!zero_TR) {
    Vmeas = Vmeas_blk.getoutput(Ec,t_inc,int_flag,true);
  } else {
    Vmeas = Ec;
  }

  double Verr = Vref - Vmeas + Vs;

  VF = Feedback_blk.getoutput(Efd,t_inc,int_flag,true);

  double Regulator_blk_in = Verr - VF;

  if(zero_TA) {
    VR = Regulator_gain_blk.getoutput(Regulator_blk_in);
  } else {
    VR = Regulator_blk.getoutput(Regulator_blk_in,t_inc,int_flag,true);
  }

  double sat_signal = KE*Efd;
  if(has_Sat) {
    sat_signal += Efd*Sat(Efd);
  }

  double output_blk_in = VR - sat_signal;

  Efd = Output_blk.getoutput(output_blk_in,t_inc,int_flag,true);
}
/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Ieeet1Model::predictor(double t_inc, bool flag)
{
  computeModel(t_inc,PREDICTOR);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Ieeet1Model::corrector(double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR);

}

/**
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void gridpack::dynamic_simulation::Ieeet1Model::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double gridpack::dynamic_simulation::Ieeet1Model::getFieldVoltage()
{
  return Efd;
}

/** 
 * Set the value of terminal voltage
 */
void gridpack::dynamic_simulation::Ieeet1Model::setVterminal(double Vm)
{
  Ec = Vm;
}

/**
 * Set stabilizer signal input
 */
void gridpack::dynamic_simulation::Ieeet1Model::setVstab(double Vstab)
{
  Vs = Vstab;
}

/** 
 * Set the exciter bus number
 */
void gridpack::dynamic_simulation::Ieeet1Model::setExtBusNum(int ExtBusNum)
{
  p_bus_num = ExtBusNum;
}	

/** 
 * Set the exciter generator id
 */
void gridpack::dynamic_simulation::Ieeet1Model::setExtGenId(std::string ExtGenId)
{
  p_gen_id = ExtGenId;
}

