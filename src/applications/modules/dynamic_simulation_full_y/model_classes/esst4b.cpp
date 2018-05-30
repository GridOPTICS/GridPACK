/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   esst4b.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   Oct 12, 2015
 * 
 * @brief  
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

#define TS_THRESHOLD 1

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Esst4bModel::Esst4bModel(void)
{
  dx1Vm = 0;
  dx2Vcomp = 0;
  dx3Va = 0; 
  dx4Vr = 0; 
  dx1Vm_1 = 0;
  dx2Vcomp_1 = 0;
  dx3Va_1 = 0;
  dx4Vr_1 = 0;
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
  //if (!data->getValue(EXCITER_KPR, &Kpr, idx)) 
  Kpr = 0.0; // TBD: Kpr
  //if (!data->getValue(EXCITER_KIR, &Kir, idx)) 
  Kir = 0.0; // TBD: Kir
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 0.0; // Vrmax
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = 0.0; // Vrmin
  if (!data->getValue(EXCITER_TA, &Ta, idx)) Ta = 0.0; // Ta
  //if (!data->getValue(EXCITER_KPM, &Kpm, idx)) 
  Kpm = 0.0; // TBD: Kpm
  //if (!data->getValue(EXCITER_KIM, &Kim, idx)) 
  Kim = 0.0; // TBD: Kim
  //if (!data->getValue(EXCITER_VMMAX, &Vmmax, idx)) 
  Vmmax = 0.0; // TBD: Vmmax
  //if (!data->getValue(EXCITER_VMMIN, &Vmmin, idx)) 
  Vmmin = 0.0; // TBD: Vmmin
  //if (!data->getValue(EXCITER_KG, &Kg, idx)) 
  Kg = 0.0; // TBD: Kg
  //if (!data->getValue(EXCITER_KP, &Kp, idx)) 
  Kp = 0.0; // TBD: Kp
  //if (!data->getValue(EXCITER_KI, &KI, idx)) 
  KI = 0.0; // TBD: KI
  //if (!data->getValue(EXCITER_VBMAX, &Vbmax, idx)) 
  Vbmax = 0.0; // TBD: Vbmax
  //if (!data->getValue(EXCITER_KC, &Kc, idx)) 
  Kc = 0.0; // TBD: Kc
  //if (!data->getValue(EXCITER_XL, &Xl, idx)) 
  Xl = 0.0; // TBD: Xl
  //if (!data->getValue(EXCITER_KPANG, &Kpang, idx)) 
  Kpang = 0.0; // TBD: Kpang
  //if (!data->getValue(EXCITER_VGMAX, &Vgmax, idx)) 
  Vgmax = 0.0; // TBD: Vgmax
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
  double result;
  if (IN <= 0.0) result = 1;
  else if (IN <= 0.433) result = 1 - 0.577 * IN;
  else if (IN <= 0.75) result = sqrt(0.75 - sqr(IN));
  else if (IN < 1.0) result = 1.732 * (1 - IN);
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
  Kpvr = Kp * cos(Kpang * 180 / pi); 
  Kpvi = Kp * sin(Kpang * 180 / pi); 
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
  Vterm = mag;
  presentMag = mag;
  Theta = ang;
  presentAng = ang;
  // State 1
  double Vb = CalculateVb(Vterm, Theta, Ir, Ii, LadIfd); // TBD: What's the init value of Ir and Ii?
  printf("esst4b: Efd = %f\n", Efd);
  double Vm = Efd/ Vb;
  // Check limits here, but these would be 
  // initial state limit violations that are not possible!
  if (OptionToModifyLimitsForInitialStateLimitViolation) { // TBD: inital value of this bool? 
    if (Vm > Vmmax) Vmmax = Vm;
    if (Vm < Vmmin) Vmmin = Vm;
  }
  double TempIn;
  if (Kim != 0) {
    x1Vm = Vm;
    TempIn = 0;
  } else {
    x1Vm = 0;
    TempIn = Vm / Kpm;
  }
  // State 3
  x3Va = Efd * Kg;
  if (x3Va > Vgmax) x3Va = Vgmax;
  x3Va = x3Va + TempIn;
  // State 4
  double Vr = x3Va;
  // Check limits here, but these would be
  // initial state limit violations that are not possible!
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (Vr > Vrmax) Vrmax = Vr;
    if (Vr < Vrmin) Vrmin = Vr;
  }
  if (Kir != 0) {
    x4Vr = Vr;
    TempIn = 0;
  } else {
    x4Vr = 0;
    TempIn = Vr / Kpr;
  }
  // State 2
  x2Vcomp = Vcomp;  // TBD: init value of Vcomp?
  // Vref
  Vref = Vcomp + TempIn;

  printf("esst4b init:  %f\t%f\t%f\t%f\n", x1Vm, x2Vcomp, x3Va, x4Vr); 
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Esst4bModel::predictor(double t_inc, bool flag)
{
  if (!flag) {
    x1Vm = x1Vm_1;
    x2Vcomp = x2Vcomp_1;
    x3Va = x3Va_1;
    x4Vr = x4Vr_1;
  }
  // State 2
  if (Tr < TS_THRESHOLD * t_inc) {
    x2Vcomp = Vcomp; // Must propogate input value instantaneously
    dx2Vcomp = 0;
  } else {
    dx2Vcomp = 1 / Tr * (Vcomp - x2Vcomp);
  }
  // State 4
  double TempIn = -x2Vcomp + Vref + Vstab;// + Vuel; // TBD: what is Vuel?
  double TempMax = Vrmax - TempIn * Kpr;
  double TempMin = Vrmin - TempIn * Kpr;
  if (x4Vr > TempMax) x4Vr = TempMax;
  else if (x4Vr < TempMin) x4Vr = TempMin;
  dx4Vr = Kir * TempIn;
  if (dx4Vr > 0 && x4Vr >= TempMax) dx4Vr = 0;
  else if (dx4Vr < 0 && x4Vr <= TempMin) dx4Vr = 0;
  //State 3
  TempIn = x4Vr + TempIn * Kpr;
  if (Ta < TS_THRESHOLD * t_inc) {
    x3Va = TempIn; // Must propogate input value instantaneously
    dx3Va = 0;
  } else {
    dx3Va = 1 / Tr * (TempIn - x3Va);
  }
  //State 1
  TempIn = Efd * Kg;
  if (TempIn > Vgmax) TempIn = Vgmax;
  TempIn = x3Va - TempIn;
  TempMax = Vmmax - TempIn * Kpm;
  TempMin = Vmmin - TempIn * Kpm;
  if (x1Vm > TempMax) x1Vm = TempMax;
  else if (x1Vm < TempMin) x1Vm = TempMin;
  dx1Vm = Kim * TempIn;
  if (dx1Vm > 0 && x1Vm >= TempMax) dx1Vm = 0;
  else if (dx1Vm < 0 && x1Vm <= TempMin) dx1Vm = 0;

  x1Vm_1 = x1Vm + dx1Vm * t_inc;
  x2Vcomp_1 = x2Vcomp + dx2Vcomp * t_inc;
  x3Va_1 = x3Va + dx3Va * t_inc;
  x4Vr_1 = x4Vr + dx4Vr * t_inc;

  printf("esst4b dx: %f\t%f\t%f\t%f\t\n", dx1Vm, dx2Vcomp, dx3Va, dx4Vr);
  printf("esst4b x: %f\t%f\t%f\t%f\n", x1Vm_1, x2Vcomp_1, x3Va_1, x4Vr_1);

  double Vb = CalculateVb(Vterm, Theta, Ir, Ii, LadIfd);
  //if (x1Vm > Voel) TempIn = Voel * Vb; // TBD: what is Voel?
  //else Efd = x1Vm_1 * Vb;
  Efd = x1Vm_1 * Vb; // TBD: temporailly

  printf("esst4b Efd: %f\n", Efd);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Esst4bModel::corrector(double t_inc, bool flag)
{
  // State 2
  if (Tr < TS_THRESHOLD * t_inc) {
    x2Vcomp_1 = Vcomp; // Must propogate input value instantaneously
    dx2Vcomp_1 = 0;
  } else {
    dx2Vcomp_1 = 1 / Tr * (Vcomp - x2Vcomp_1);
  }
  // State 4
  double TempIn = -x2Vcomp_1 + Vref + Vstab;// + Vuel; // TBD: what is Vuel?
  double TempMax = Vrmax - TempIn * Kpr;
  double TempMin = Vrmin - TempIn * Kpr;
  if (x4Vr_1 > TempMax) x4Vr_1 = TempMax;
  else if (x4Vr_1 < TempMin) x4Vr_1 = TempMin;
  dx4Vr_1 = Kir * TempIn;
  if (dx4Vr_1 > 0 && x4Vr_1 >= TempMax) dx4Vr_1 = 0;
  else if (dx4Vr_1 < 0 && x4Vr_1 <= TempMin) dx4Vr_1 = 0;
  //State 3
  TempIn = x4Vr_1 + TempIn * Kpr;
  if (Ta < TS_THRESHOLD * t_inc) {
    x3Va_1 = TempIn; // Must propogate input value instantaneously
    dx3Va_1 = 0;
  } else {
    dx3Va_1 = 1 / Tr * (TempIn - x3Va_1);
  }
  //State 1
  TempIn = Efd * Kg;
  if (TempIn > Vgmax) TempIn = Vgmax;
  TempIn = x3Va_1 - TempIn;
  TempMax = Vmmax - TempIn * Kpm;
  TempMin = Vmmin - TempIn * Kpm;
  if (x1Vm_1 > TempMax) x1Vm_1 = TempMax;
  else if (x1Vm_1 < TempMin) x1Vm_1 = TempMin;
  dx1Vm_1 = Kim * TempIn;
  if (dx1Vm_1 > 0 && x1Vm_1 >= TempMax) dx1Vm_1 = 0;
  else if (dx1Vm_1 < 0 && x1Vm_1 <= TempMin) dx1Vm_1 = 0;

  x1Vm_1 = x1Vm + (dx1Vm + dx1Vm_1) / 2.0 * t_inc;
  x2Vcomp_1 = x2Vcomp + (dx2Vcomp + dx2Vcomp_1) / 2.0 * t_inc;
  x3Va_1 = x3Va + (dx3Va + dx3Va_1) / 2.0 * t_inc;
  x4Vr_1 = x4Vr + (dx4Vr + dx4Vr_1) / 2.0 * t_inc;

  printf("esst4b dx: %f\t%f\t%f\t%f\t\n", dx1Vm_1, dx2Vcomp_1, dx3Va_1, dx4Vr_1);
  printf("esst4b x: %f\t%f\t%f\t%f\n", x1Vm_1, x2Vcomp_1, x3Va_1, x4Vr_1);
  
  double Vb = CalculateVb(Vterm, Theta, Ir, Ii, LadIfd);
  //if (x1Vm > Voel) TempIn = Voel * Vb; // TBD: what is Voel?
  //else Efd = x1Vm_1 * Vb;
  Efd = x1Vm_1 * Vb; // TBD: temporially

  printf("esst4b Efd: %f\n", Efd);
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
  //Vterminal = mag;
}

/** 
 * Set the value of the omega
 * @return value of field current
 */
void gridpack::dynamic_simulation::Esst4bModel::setOmega(double omega)
{
  //w = omega;
}

