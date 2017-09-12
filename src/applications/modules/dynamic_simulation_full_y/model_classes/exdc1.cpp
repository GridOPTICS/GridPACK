/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   exdc1.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 11, 2015
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_exciter_model.hpp"
#include "exdc1.hpp"

#define TS_THRESHOLD 1

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Exdc1Model::Exdc1Model(void)
{
  dx1 = 0;
  dx2 = 0;
  dx3 = 0; 
  dx4 = 0; 
  dx5 = 0;
  dx1_1 = 0;
  dx2_1 = 0;
  dx3_1 = 0;
  dx4_1 = 0;
  dx5_1 = 0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Exdc1Model::~Exdc1Model(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * Exdc1Model
 */
void gridpack::dynamic_simulation::Exdc1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(EXCITER_TR, &TR, idx)) TR = 0.0; // TR
  if (!data->getValue(EXCITER_KA, &KA, idx)) KA = 0.0; // KA 
  if (!data->getValue(EXCITER_TA, &TA, idx)) TA = 0.0; // TA
  if (!data->getValue(EXCITER_TB, &TB, idx)) TB = 0.0; // TB
  if (!data->getValue(EXCITER_TC, &TC, idx)) TC = 0.0; // TC
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 0.0; // Vrmax
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = 0.0; // Vrmin
  if (!data->getValue(EXCITER_KE, &KE, idx)) KE = 0.0; // KE
  if (!data->getValue(EXCITER_TE, &TE, idx)) TE = 0.0; // TE
  if (!data->getValue(EXCITER_KF, &KF, idx)) KF = 0.0; // KF
  //if (!data->getValue(EXCITER_TF1, &TF1, idx)) TF1 = 0.0; // TF1
  if (!data->getValue(EXCITER_TF1, &TF, idx)) TF = 0.0; // TF
  //printf("load TF = %f\n", TF);
  if (!data->getValue(EXCITER_SWITCH, &SWITCH, idx)) SWITCH = 0.0; // SWITCH
  if (!data->getValue(EXCITER_E1, &E1, idx)) E1 = 0.0; // E1
  if (!data->getValue(EXCITER_SE1, &SE1, idx)) SE1 = 0.0; // SE1
  if (!data->getValue(EXCITER_E2, &E2, idx)) E2 = 0.0; // E2
  if (!data->getValue(EXCITER_SE2, &SE2, idx)) SE2 = 0.0; // SE2
  //if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.0; // S10
  //if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.0; // S12
}

/**
 *  * Saturation function
 *   * @ param x
 *    */
double gridpack::dynamic_simulation::Exdc1Model::Sat(double x)
{
    /*double a_ = S12 / S10 - 1.0 / 1.2;
    double b_ = (-2 * S12 / S10 + 2);
    double c_ = S12 / S10 - 1.2;
    double A = (-b_ - sqrt(b_ * b_ - 4 * a_ * c_)) / (2 * a_);
    double B = S10 / ((1.0 - A) * (1.0 - A));
    return B * ( x - A) * (x - A);*/
    double B = log(SE2 / SE1)/(E2 - E1);
    double A = SE1 / exp(B * E1);
    return A * exp(B * x);
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Exdc1Model::init(double mag, double ang, double ts)
{
  ///printf("exdc1: Efd = %f\n", Efd);
  x1 = Efd;
  x4 = Efd * (KE + Sat(Efd));
  if (TB > (TS_THRESHOLD * ts)) 
    x3 = (x4 / KA) * (1 - TC / TB); // SJin: x4 is Vr 
  else
    x3 = x4 / KA;
  x2 = mag; // SJin: mag is Vterminal 
  //printf("KF = %f, TF = %f, x1 = %f, ts = %f\n", KF, TF, x1, ts);
  if (TF > (TS_THRESHOLD * ts)) 
    x5 = x1 * (KF / TF); // SJin: x1 is Ve
  else
    x5 = 0.0;
  //x5 = 0.0; // Diao: force x5 to 0.0 for now
  Vref = mag + x4 / KA;
  //printf("Vref = %f\n", Vref);
  ///printf("exdc1 init:  %f\t%f\t%f\t%f\t%f\n", x1, x2, x3, x4, x5); 
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Exdc1Model::predictor(double t_inc, bool flag)
{
  if (!flag) {
    x1 = x1_1;
    x2 = x2_1;
    x3 = x3_1;
    x4 = x4_1;
    x5 = x5_1;
  }
  //printf("...........%f\t%f\t%f\t%f\t%f\n", x1, x2, x3, x4, x5);
  //printf(".........Vterminal = %f\n", Vterminal);
  double Feedback;
  // State 2
  //printf("TR = %f\n", TR);
  if (TR > TS_THRESHOLD * t_inc) dx2 = (Vterminal - x2) / TR; // SJin: x2 is Vsens
  else x2 = Vterminal; // propogate state immediately
  // State 5
  if (TF > TS_THRESHOLD * t_inc) {
    dx5 = (x1 * KF / TF - x5) / TF;
    Feedback = x1 * KF / TF - x5;
  } else {
    x5 = 0;
    Feedback = 0;
  }
  //printf("TF = %f, t_inc = %f, x5 = %f, dx5 = %f\n", TF, t_inc, x5, dx5);
  // State 3
  double Vstab = 0.0; // SJin: Output from PSS, set to 0.0 for now.
  double LeadLagIN = Vref - x2 + Vstab - Feedback;
  //printf("Vref = %f, TB = %f, TC = %f\n", Vref, TB, TC); 
  double LeadLagOUT;
  //printf("LeadLagIN = %f, TC = %f, TB = %f, x3 = %f\n", LeadLagIN, TC, TB, x3);
  if (TB > (TS_THRESHOLD * t_inc)) {
    dx3 = (LeadLagIN * (1 - TC / TB) - x3) / TB;
    LeadLagOUT = LeadLagIN * TC / TB + x3;
  } else 
    LeadLagOUT = LeadLagIN;
  // State 4
  //printf("x4 = %f, Vrmax = %f, Vrmin = %f, TA = %f, LeadLagOUT = %f, KA = %f\n", x4, Vrmax, Vrmin, TA, LeadLagOUT, KA);
  if (x4 > Vrmax) x4 = Vrmax;
  if (x4 < Vrmin) x4 = Vrmin;
  if (TA > (TS_THRESHOLD * t_inc)) 
    dx4 = (LeadLagOUT * KA - x4) / TA;
  else {
    dx4 = 0;
    x4 = LeadLagOUT * KA;
  }
  if (dx4 > 0 && x4 >= Vrmax) dx4 = 0;
  if (dx4 < 0 && x4 <= Vrmin) dx4 = 0;
  // State 1
  dx1 = x4 - x1 * (KE + Sat(x1));

  x1_1 = x1 + dx1 * t_inc;
  x2_1 = x2 + dx2 * t_inc;
  x3_1 = x3 + dx3 * t_inc;
  x4_1 = x4 + dx4 * t_inc;
  x5_1 = x5 + dx5 * t_inc;
  //printf("x2 = %f, dx2 = %f, x2_1 = %f\n", x2, dx2, x2_1);
  //printf("x5 = %f, dx5 = %f, x5_1 = %f\n", x5, dx5, x5_1);

  ///printf("exdc1 dx: %f\t%f\t%f\t%f\t%f\t\n", dx1, dx2, dx3, dx4, dx5);
  ///printf("exdc1 x: %f\t%f\t%f\t%f\t%f\n", x1_1, x2_1, x3_1, x4_1, x5_1);

  Efd = x1_1 * (1 + w);

  ///printf("exdc1 Efd: %f\n", Efd);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Exdc1Model::corrector(double t_inc, bool flag)
{
  double Feedback;
  // State 2
  //printf("TR = %f, Vterminal = %f, x2_1 = %f\n", TR, Vterminal, x2_1);
  if (TR > TS_THRESHOLD * t_inc) dx2_1 = (Vterminal - x2_1) / TR; // SJin: x2 is Vsens
  else x2_1 = Vterminal; // propogate state immediately
  // State 5
  if (TF > TS_THRESHOLD * t_inc) {
    dx5_1 = (x1_1 * KF / TF - x5_1) / TF;
    Feedback = x1_1 * KF / TF - x5_1;
  } else {
    x5_1 = 0;
    Feedback = 0;
  }
  // State 3
  double Vstab = 0.0; // SJin: Output from PSS, set to 0.0 for now.
  double LeadLagIN = Vref - x2_1 + Vstab - Feedback;  
  double LeadLagOUT;
  //printf("LeadLagIN = %f, TC = %f, TB = %f, x3 = %f\n", LeadLagIN, TC, TB, x3);
  if (TB > (TS_THRESHOLD * t_inc)) {
    dx3_1 = (LeadLagIN * (1 - TC / TB) - x3_1) / TB;
    LeadLagOUT = LeadLagIN * TC / TB + x3_1;
  } else 
    LeadLagOUT = LeadLagIN;
  // State 4
  if (x4_1 > Vrmax) x4_1 = Vrmax;
  if (x4_1 < Vrmin) x4_1 = Vrmin;
  if (TA > (TS_THRESHOLD * t_inc)) 
    dx4_1 = (LeadLagOUT * KA - x4_1) / TA;
  else {
    dx4_1 = 0;
    x4_1 = LeadLagOUT * KA;
  }
  if (dx4_1 > 0 && x4_1 >= Vrmax) dx4_1 = 0;
  if (dx4_1 < 0 && x4_1 <= Vrmin) dx4_1 = 0;
  // State 1
  dx1_1 = x4_1 - x1_1 * (KE + Sat(x1));

  x1_1 = x1 + (dx1 + dx1_1) / 2.0 * t_inc;
  x2_1 = x2 + (dx2 + dx2_1) / 2.0 * t_inc;
  x3_1 = x3 + (dx3 + dx3_1) / 2.0 * t_inc;
  x4_1 = x4 + (dx4 + dx4_1) / 2.0 * t_inc;
  x5_1 = x5 + (dx5 + dx5_1) / 2.0 * t_inc;

  ///printf("exdc1 dx: %f\t%f\t%f\t%f\t%f\t\n", dx1_1, dx2_1, dx3_1, dx4_1, dx5_1);
  ///printf("exdc1 x: %f\t%f\t%f\t%f\t%f\n", x1_1, x2_1, x3_1, x4_1, x5_1);
  
  Efd = x1_1 * (1 + w);

  ///printf("exdc1 Efd: %f\n", Efd);
}

/**
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void gridpack::dynamic_simulation::Exdc1Model::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

/**
 * Set the field current parameter inside the exciter
 * @param fldc value of the field current
 */
void gridpack::dynamic_simulation::Exdc1Model::setFieldCurrent(double fldc)
{
  LadIfd = fldc;
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double gridpack::dynamic_simulation::Exdc1Model::getFieldVoltage()
{
  return Efd;
}

/** 
 * Get the value of the field current parameter
 * @return value of field current
 */
double gridpack::dynamic_simulation::Exdc1Model::getFieldCurrent()
{
  return 0.0;
}

/** 
 * Set the value of the Vterminal
 * @return value of field current
 */
void gridpack::dynamic_simulation::Exdc1Model::setVterminal(double mag)
{
  Vterminal = mag;
}

/** 
 * Set the value of the omega
 * @return value of field current
 */
void gridpack::dynamic_simulation::Exdc1Model::setOmega(double omega)
{
  w = omega;
}

