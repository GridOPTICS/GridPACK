/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   esst4b.cpp
 * 
 * @brief  ESST4B
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_exciter_model.hpp"
#include "esst4b.hpp"

#define TS_THRESHOLD 4
#define KI_THRESHOLD 1e-6

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Esst4bModel::Esst4bModel(void)
{

}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Esst4bModel::~Esst4bModel(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * Esst4bModel
 */
void gridpack::dynamic_simulation::Esst4bModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(EXCITER_TR, &Tr, idx)) Tr = 0.0; // Tr
  if (!data->getValue(EXCITER_KPR, &Kpr, idx)) Kpr = 0.0; // Kpr
  if (!data->getValue(EXCITER_KIR, &Kir, idx)) Kir = 0.0; // Kir
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 0.0; // Vrmax
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = 0.0; // Vrmin
  if (!data->getValue(EXCITER_TA, &Ta, idx)) Ta = 0.0; // Ta
  if (!data->getValue(EXCITER_KPM, &Kpm, idx)) Kpm = 0.0; // Kpm
  if (!data->getValue(EXCITER_KPM, &Kim, idx)) Kim = 0.0; // Kim
  if (!data->getValue(EXCITER_VMMAX, &Vmmax, idx)) Vmmax = 0.0; // Vmmax
  if (!data->getValue(EXCITER_VMMIN, &Vmmin, idx)) Vmmin = 0.0; // Vmax
  if (!data->getValue(EXCITER_KG, &Kg, idx)) Kg = 0.0; // Kg
  if (!data->getValue(EXCITER_KP, &Kp, idx)) Kp = 0.0; // Kp
  if (!data->getValue(EXCITER_KI, &KI, idx)) KI = 0.0; // KI
  if (!data->getValue(EXCITER_VBMAX, &Vbmax, idx)) Vbmax = 0.0; // Vbmax
  if (!data->getValue(EXCITER_KC, &Kc, idx)) Kc = 0.0; // Kc
  if (!data->getValue(EXCITER_XL, &Xl, idx)) Xl = 0.0; // Xl
  if (!data->getValue(EXCITER_THETAP, &Thetap, idx)) Thetap = 0.0; // Thetap, yuan: deg
  

  // right now we just hard code Vuel, Voel and Vs
  Vs = 0.0;
  Vuel = 0.0;
  Voel = 1000.0;
  OptionToModifyLimitsForInitialStateLimitViolation = true;
  zero_TR = false;
  zero_TA = false;
  zero_KIM = false;
  zero_KIR = false;

}

/**
 * Saturation function
 * @ param x
 */
double gridpack::dynamic_simulation::Esst4bModel::Sat(double x)
{
    /*double a_ = S12 / S10 - 1.0 / 1.2;
    double b_ = (-2 * S12 / S10 + 2);
    double c_ = S12 / S10 - 1.2;
    double A = (-b_ - sqrt(b_ * b_ - 4 * a_ * c_)) / (2 * a_);
    double B = S10 / ((1.0 - A) * (1.0 - A));
    return B * ( x - A) * (x - A);*/
    //double B = log(SE2 / SE1)/(E2 - E1);
    //double A = SE1 / exp(B * E1);
    //return A * exp(B * x);

    return 0.0;
}

double gridpack::dynamic_simulation::Esst4bModel::sqr(double x)
{
  return x * x;
}

/**
 * FEX function
 * @ param IN
 */
double gridpack::dynamic_simulation::Esst4bModel::FEX(double IN)
{

  // double result;
  // if (IN <= 0.0) result = 1;
  // else if (IN <= 0.433) result = 1 - 0.577 * IN;
  // else if (IN <= 0.75) result = sqrt(0.75 - sqr(IN));
  // else if (IN < 1.0) result = 1.732 * (1 - IN);
  // else result = 0;
  // return result;
  
  // the above old code might be PowerWorld implementation; the new was referred to PSS/E 35 model library
  // still different from PSS/E implementation, a little confusing here.
  double result;
  if (IN <= 0.0) result = 1;
  else if (IN <= 0.433 && IN > 0.0) result = 1 - 0.577 * IN;
  else if (IN < 0.75 && IN > 0.433) result = sqrt(0.75 - sqr(IN));
  else if (IN <= 1.0 && IN >= 0.75) result = 1.732 * (1 - IN);
  else result = 0;
  return result;
  
}

/**
 * CalculateVb function
 * @ param Vterm, Theta, Ir, Ii, LadIfd
 */
double gridpack::dynamic_simulation::Esst4bModel::CalculateVb(double Vterm,
  double Theta, double Ir, double Ii, double LadIfd)
{
  double pi = 4.0 * atan(1.0);
  // Calculate complex values for VE Calculation
  Kpvr = Kp * cos(Thetap * pi / 180.0); // Thetap and Kp were read in from Parser in the original implementation.
  Kpvi = Kp * sin(Thetap * pi / 180.0); // cos and sin here take rad as inputs; Thetap is in deg
  Kpir = - Kpvi * Xl;
  Kpii = + Kpvr * Xl + KI;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Ve = sqr(Kpvr * Vrterm - Kpvi * Viterm + Kpir * Ir - Kpii * Ii) +
              sqr(Kpvr * Viterm + Kpvi * Vrterm + Kpir * Ii + Kpii * Ir);
  Ve = sqrt(Ve);
  double Vb = Ve * FEX(Kc * LadIfd / Ve); // FEX function from 1.8.1
  if (Vb > Vbmax) Vb = Vbmax;
  return Vb; 
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Esst4bModel::init(double mag, double ang, double ts)
{
  /*---yuan add below---*/
  if (Tr < TS_THRESHOLD * ts) zero_TR = true;
  if (Ta < TS_THRESHOLD * ts) zero_TA = true;
  if (Kim < KI_THRESHOLD) zero_KIM = true;
  if (Kir < KI_THRESHOLD) zero_KIR = true;
  
  if (zero_TR) printf("Tr=%f is better at least %d times larger than timestep=%f.\n", Tr, TS_THRESHOLD, ts);
  if (zero_TA) printf("Ta=%f is better at least %d times larger than timestep=%f.\n", Ta, TS_THRESHOLD, ts);
  if (zero_KIM) printf("Kim=%f is less than %d times of timestep=%f, treated as zero.\n", Kim, TS_THRESHOLD, ts);
  if (zero_KIR) printf("Kir=%f is less than %d times of timestep=%f, treated as zero.\n", Kir, TS_THRESHOLD, ts);
  
  // printf("print: inside esst4b model, Ir=%f, Ii=%f\n", Ir, Ii);
  
  if (Kpr == 0.0 && Kir < KI_THRESHOLD) {
    Kpr = 1.0;
    printf("Kpr and Kir cannot be both zeros; reset Kpr=1.0.\n");
  }
  
  if (Kpm == 0.0 && Kim < KI_THRESHOLD) {
    Kpm = 1.0;
    printf("Kpm and Kim cannot be both zeros; reset Kpm=1.0.\n");
  }

  if (!zero_TR) {
    Filter_blkR.setparams(1.0, Tr);
  }
  /*---else {a gain=1 block}---*/

  if (!zero_KIR) {
    PIControl_blkR.setparams(Kpr, Kir, Vrmin, Vrmax, -10000.0, 10000.0);
  }
  /*---else {a gain=Kpr block}---*/
  
  if (!zero_TA) {
    Filter_blkA.setparams(1.0, Ta);
  }
  
  if (!zero_KIM) {
    PIControl_blmM.setparams(Kpm, Kim, Vmmin, Vmmax, -10000.0, 10000.0);
  }
  
  LVGate_blk.setparams(Voel);
  
  double u1, u2, u3, u4;

  Vterm = mag; // Ec
  Theta = ang;

  double Vb = CalculateVb(Vterm, Theta, Ir, Ii, LadIfd); 

  u1 = Efd / Vb;

  // LV Gate?
  // assume LV gate always take feedforward path during initialization, may need to adjust Voel
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (u1 > Voel) Voel = u1 + 0.1;
    LVGate_blk.setparams(Voel);
  }
  
  
  if (!zero_KIM) {
    u2 = PIControl_blmM.init_given_y(u1);
    // Check limits here, but these would be 
    // initial state limit violations that are not possible!
    if (OptionToModifyLimitsForInitialStateLimitViolation) {
      if (u1 > Vmmax) Vmmax = u1+0.1;
      if (u1 < Vmmin) Vmmin = u1-0.1;
    }
  } else {
    u2 = u1/Kpm;
  }
  
  if (!zero_TA) {
    u3 = Filter_blkA.init_given_y(u2 + Efd * Kg);
  } else {
    u3 = u2 + Efd * Kg;
  }
  
  if (!zero_KIR) {
    u4 = PIControl_blkR.init_given_y(u3);
    // Check limits here, but these would be 
    // initial state limit violations that are not possible!
    if (OptionToModifyLimitsForInitialStateLimitViolation) {
      if (u3 > Vrmax) Vrmax = u3+0.1;
      if (u3 < Vrmin) Vrmin = u3-0.1;
    }
  } else {
    u4 = u3/Kpr;
  }

  if(!zero_TR) {
    Vmeas = Filter_blkR.init_given_u(Vcomp);
  } else 
    Vmeas = Vcomp;
  Vref = u4 + Vmeas - Vuel - Vs; 

}

/**
 * computeModel - Updates the model states and gets the output
 */
void gridpack::dynamic_simulation::Esst4bModel::computeModel(double t_inc,IntegrationStage int_flag)
{
  double u1, u2, u3, u4;
  
  
  if(!zero_TR) {
    Vmeas = Filter_blkR.getoutput(Vcomp, t_inc, int_flag, true);
  } else {
    Vmeas = Vcomp;
  }
  u1 = Vref - Vmeas + Vuel + Vs;
  
  if(!zero_KIR) {
    u2 = PIControl_blkR.getoutput(u1, t_inc, int_flag, true);
  } else {
    u2 = Kpr * u1;
  }
  
  if(!zero_TA) {
    u3 = Filter_blkA.getoutput(u2, t_inc, int_flag, true);
  } else {
    u3 = u2;
  }
  
  if(!zero_KIM) {
    u4 = PIControl_blmM.getoutput(u3 - Efd * Kg, t_inc, int_flag, true);
  } else {
    u4 = Kpm * (u3 - Efd * Kg);
  }
  u4 = LVGate_blk.getoutput(u4);
  double Vb = CalculateVb(Vterm, Theta, Ir, Ii, LadIfd); 
  Efd = u4 * Vb;

}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Esst4bModel::predictor(double t_inc, bool flag)
{
  computeModel(t_inc,PREDICTOR);
  
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Esst4bModel::corrector(double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR);

}

/**
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void gridpack::dynamic_simulation::Esst4bModel::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

/**
 * Set the field current parameter inside the exciter
 * @param fldc value of the field current
 */
void gridpack::dynamic_simulation::Esst4bModel::setFieldCurrent(double fldc)
{
  LadIfd = fldc;
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double gridpack::dynamic_simulation::Esst4bModel::getFieldVoltage()
{
  return Efd;
}

/** 
 * Get the value of the field current parameter
 * @return value of field current
 */
double gridpack::dynamic_simulation::Esst4bModel::getFieldCurrent()
{
  return 0.0;
}

/** 
 * Set the value of the Vterminal
 * @return value of field current
 */
void gridpack::dynamic_simulation::Esst4bModel::setVterminal(double mag)
{
  Vterm = mag;
}

/** 
 * Set the value of the omega
 * @return value of field current
 */
void gridpack::dynamic_simulation::Esst4bModel::setOmega(double omega)
{
  //w = omega;
}

void gridpack::dynamic_simulation::Esst4bModel::setVuel(double vtmp)
{
  Vuel = vtmp;
}

void gridpack::dynamic_simulation::Esst4bModel::setVs(double vtmp)
{
  Vs = vtmp;
}

void gridpack::dynamic_simulation::Esst4bModel::setIri(double vIr, double vIi)
{
  Ir = vIr;
  Ii = vIi;
}

void gridpack::dynamic_simulation::Esst4bModel::setVoel(double vtmp)
{
  Voel = vtmp;
}

/** 
 * Set the value of the Vcomp
 * @return value of the Vcomp
 */
void gridpack::dynamic_simulation::Esst4bModel::setVcomp(double vtmp)
{
  Vcomp = vtmp;
}
