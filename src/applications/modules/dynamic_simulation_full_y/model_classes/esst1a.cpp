/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   esst1a.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   Oct 12, 2015
 * @Last modified with control block: Sep 1, 2023
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
#include "esst1a.hpp"

#define TS_THRESHOLD 1

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Esst1aModel::Esst1aModel(void)
{
  /*dx1Va = 0;
  dx2Vcomp = 0;
  dx3LL1 = 0; 
  dx4LL2 = 0; 
  dx5Deriv = 0;
  dx1Va_1 = 0;
  dx2Vcomp_1 = 0;
  dx3LL1_1 = 0;
  dx4LL2_1 = 0;
  dx5Deriv_1 = 0;
  OptionToModifyLimitsForInitialStateLimitViolation = true;*/

  zero_TA = false;
  zero_TR = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Esst1aModel::~Esst1aModel(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * Esst1aModel
 */
void gridpack::dynamic_simulation::Esst1aModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(EXCITER_UEL, &UEL, idx)) UEL = 0.0; // TBD: UEL
  if (!data->getValue(EXCITER_VOS, &VOS, idx)) VOS = 0.0; // TBD: VOS
  if (!data->getValue(EXCITER_TR, &Tr, idx)) Tr = 0.0; // Tr
  if (!data->getValue(EXCITER_VIMAX, &Vimax, idx)) Vimax = 0.0; // TBD: Vimax
  if (!data->getValue(EXCITER_VIMIN, &Vimin, idx)) Vimin = 0.0; // TBD: Vimin
  if (!data->getValue(EXCITER_TC, &Tc, idx)) Tc = 0.0; // Tc
  if (!data->getValue(EXCITER_TB, &Tb, idx)) Tb = 0.0; // Tb
  if (!data->getValue(EXCITER_TC1, &Tc1, idx)) Tc1 = 0.0; // TBD: Tc1
  if (!data->getValue(EXCITER_TB1, &Tb1, idx)) Tb1 = 0.0; // TBD: Tb1
  if (!data->getValue(EXCITER_KA, &Ka, idx)) Ka = 0.0; // Ka
  if (!data->getValue(EXCITER_TA, &Ta, idx)) Ta = 0.0; // Ta
  if (!data->getValue(EXCITER_VAMAX, &Vamax, idx)) Vamax = 0.0; // TBD: Vamax
  if (!data->getValue(EXCITER_VAMIN, &Vamin, idx)) Vamin = 0.0; // TBD: Vamin
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 0.0; // Vrmax
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = 0.0; // Vrmin
  if (!data->getValue(EXCITER_KC, &Kc, idx)) Kc = 0.0; // TBD: Kc
  if (!data->getValue(EXCITER_KF, &Kf, idx)) Kf = 0.0; // Kf
  if (!data->getValue(EXCITER_TF, &Tf, idx)) Tf = 0.0; // Tf
  if (!data->getValue(EXCITER_KLR, &Klr, idx)) Klr = 0.0; // TBD: Klr
  if (!data->getValue(EXCITER_ILR, &Ilr, idx)) Ilr = 0.0; // TBD: Ilr

  //printf("UEL=%f,VOS=%f,Tr=%f,Vimax=%f,Vimin=%f,Tc=%f,Tb=%f,Tc1=%f,Tb1=%f,Ka=%f,Ta=%f,Vamax=%f,Vamin=%f,Vrmax=%f,Vrmin=%f,Kc=%f,Kf=%f,Tf=%f,Klr=%f,Ilr=%f\n",UEL,VOS,Tr,Vimax,Vimin,Tc,Tb,Tc1,Tb1,Ka,Ta,Vamax,Vamin,Vrmax,Vrmin,Kc,Kf,Tf,Klr,Ilr);

  if(fabs(Ta) < 1e-6) zero_TA = true;
  if(fabs(Tr) < 1e-6) zero_TR = true;

  if(!zero_TR) {
    Filter_blkR.setparams(1.0, Tr);
  }

  HVGate_blk1.setparams(Vuel); // is UEL Vuel?

  Leadlag_blkBC.setparams(Tc, Tb);
  Leadlag_blkBC1.setparams(Tc1, Tb1);

  if(!zero_TA) {
    Regulator_blk.setparams(Ka,Ta,Vrmin,Vrmax,-1000.0,1000.0);
  } else {
    Regulator_gain_blk.setparams(Ka,Vrmin,Vrmax);
  }

  HVGate_blk2.setparams(Vuel); // UEL is Vuel?
  LVGate_blk.setparams(Voel); // Where to read Voel from? Set it by funciton call as Vuel?

  double a[2], b[2];
  a[0] = Tf; a[1] = 1.0;
  b[0] = Kf; b[1] = 0.0;
  Feedback_blkF.setcoeffs(a, b);

}

/**
 * Saturation function
 * @ param x
 */
double gridpack::dynamic_simulation::Esst1aModel::Sat(double x)
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

double gridpack::dynamic_simulation::Esst1aModel::sqr(double x)
{
  return x * x;
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Esst1aModel::init(double mag, double ang, double ts)
{
  Vf = Feedback_blkF.init_given_u(Efd);
  
  // LV Gate?

  // HV Gate?

  VA = (LadIfd - Ilr) * Klr - Vothsg - Efd;
  
  if (zero_TA) {
    VLL1 = VA/Ka;
  } else {
    VLL1 = Regulator_blk.init_given_y(VA);
  }

  VLL = Leadlag_blkBC1.init_given_y(VLL1);
  double u1 = Leadlag_blkBC.init_given_y(VLL);

  // HV Gate?

  double Verr = u1 - Vf + Vothsg + Vuel; 

  if(!zero_TR) {
    Vmeas = Filter_blkR.init_given_u(Vterm);
  } else 
    Vmeas = Vterm;

  Vref = Verr + Vmeas; 

 /* // Parameter cleanup
  // Just to make the code simpler below we will do the following cleanup 
  // on the input parameters and also store some variables.
  if (Tb == 0) Tc = 0;
  if (Tb1 == 0) Tc1 = 0;
  if (Tf == 0) Kf = 0;
  // Following is needed to avoid an algebraic loop in the model
  if (Kf != 0 && Ta < TS_THRESHOLD * ts) Ta = TS_THRESHOLD * ts;

  Vterm = mag;
  presentMag = mag;

  //State 1
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (Efd > (Vterm * Vrmax - Kc * LadIfd)) Vrmax = (Efd + Kc * LadIfd) / Vterm+0.21;
    if (Efd < (Vterm * Vrmin)) Vrmin = Efd / Vterm-0.1;
  }
  double Va;
  if (LadIfd > Ilr) Va = Efd + Klr * (LadIfd - Ilr);
  else Va = Efd;
  // Check for initial limit violations - should not happen!
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (Va > Vamax) Vamax = Va+0.1;
    if (Va < Vamin) Vamin = Va-0.1;
  }
  x1Va = Va;
  //State 4
  double TempLL = Va / Ka;
  if (Tb1 == 0) x4LL2 = TempLL;
  else if (Tb1 < TS_THRESHOLD * ts) x4LL2 = 0;
  else x4LL2 = TempLL * (1 - Tc1 / Tb1);
  //State 3
  if (Tb == 0) x3LL1 = TempLL;
  else if (Tb < TS_THRESHOLD * ts) x3LL1 = 0;
  else x3LL1 = TempLL * (1 - Tc / Tb);
  //State 2
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (TempLL > Vimax) Vimax = TempLL+0.1;
    if (TempLL < Vimin) Vimin = TempLL-0.1;
  }
  printf ("Esst1a ini() Vcomp = %f, TempLL = %f \n", Vcomp, TempLL);
  x2Vcomp = Vcomp;
  //State 5
  x5Deriv = 0;
  Vref = Vcomp + TempLL;
  Vstab = 0.0;

  //printf("test: esst1a init:  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f \n", x1Va, x2Vcomp, x3LL1, x4LL2, x5Deriv); 
  */
}

/**
 * computeModel - Updates the model states and gets the output
 */
void gridpack::dynamic_simulation::Esst1aModel::computeModel(double t_inc,IntegrationStage int_flag)
{
  if(!zero_TR) {
    Vmeas = Filter_blkR.getoutput(Vterm, t_inc, int_flag, true);
  } else {
    Vmeas = Vterm;
  }

  double Verr = Vref - Vmeas + Vothsg + Vuel; 

  Vf = Feedback_blkF.getoutput(Efd, t_inc, int_flag, true);

  double leadlag_blk_in = Verr - Vf;

  if (leadlag_blk_in > Vimax)
    leadlag_blk_in = Vimax;
  else if (leadlag_blk_in < Vimin)
    leadlag_blk_in = Vimin;

  leadlag_blk_in = HVGate_blk1.getoutput(leadlag_blk_in);

  VLL = Leadlag_blkBC.getoutput(leadlag_blk_in, t_inc, int_flag, true);
  
  VLL1 = Leadlag_blkBC1.getoutput(VLL, t_inc, int_flag, true); 

  if (zero_TA) {
    VA = Regulator_gain_blk.getoutput(VLL1);
  } else {
    VA = Regulator_blk.getoutput(VLL, t_inc, int_flag, true);
  }

  double u1 = VA + Vothsg - (LadIfd - Ilr) * Klr; 

  double u2 = HVGate_blk2.getoutput(u1);

  u2 = LVGate_blk.getoutput(u2);

  double VT; // to be updated since VT is not a parameter that has been given.
  if (u2 > VT * Vrmax - Kc * LadIfd)
    u2 = VT * Vrmax - Kc * LadIfd;
  else if (u2 < VT * Vrmin)
    u2 = VT * Vrmin;

  Efd = u2; 

}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Esst1aModel::predictor(double t_inc, bool flag)
{
  computeModel(t_inc,PREDICTOR);
 /* //printf ("esst1a predictor Vterm=%f, Vcomp=%f, LadIfd=%f\n", Vterm, Vcomp, LadIfd);
  if (!flag) {
    x1Va = x1Va_1;
    x2Vcomp = x2Vcomp_1;
    x3LL1 = x3LL1_1;
    x4LL2 = x4LL2_1;
    x5Deriv = x5Deriv_1;
  }
  printf("\n esst1a: what's the initial values for the first iteration?\n");
  printf(" %f\t%f\t%f\t%f\t%f\t%f\n", x1Va, x2Vcomp, x3LL1, x4LL2, x5Deriv, Vref);
  // State 2
  if (Tr < TS_THRESHOLD * t_inc) {
    x2Vcomp = Vcomp; // Must propogate input value instantaneously
    dx2Vcomp = 0;
  } else {
    dx2Vcomp = 1 / Tr * (Vcomp - x2Vcomp);
  }
  // State 5
  double TempIn;
  if (Kf > 0) {
    //TempIn = Klr * (LadIfd - Ilr);
    //if (TempIn < 0) TempIn = 0;
    // Note: if Ta = 0, then we would have an algebriac loop here which
    // could cause numerical troubles. Thus we enforce Ta > 0 if (Kf > 0 and Tf > 0)
    //TempIn = x1Va - TempIn;
    //if ("VOS at Output") TempIn = TempIn + Vstab; // TBD: "VOS at Output"???
    // Ignore Over and Under Excitation Limit for now
	TempIn = x1Va;
	
    double UseTf;
    if (Tf > TS_THRESHOLD * t_inc) UseTf = Tf;
    else UseTf = TS_THRESHOLD * t_inc;
    dx5Deriv = (TempIn * (-Kf / UseTf) - x5Deriv) / UseTf;
    TempIn = TempIn * Kf / UseTf + x5Deriv; // Input into summation block
  } else {
    dx5Deriv = 0;
    TempIn = 0; // Input into summation block;
  }
  // State 3
  TempIn = - x2Vcomp - TempIn + Vref;
  TempIn = TempIn + Vstab;
  //if ("VOS at Input") TempIn = TempIn + Vstab; // TBD: "VOS at Input"???
  if (TempIn > Vimax) TempIn = Vimax;
  if (TempIn < Vimin) TempIn = Vimin;
  // Ignore Under Excitation Limit
  if (Tb < TS_THRESHOLD * t_inc) {
    dx3LL1 = 0;
    TempIn = TempIn; // output of first lead lag - just pass input
  } else {
    dx3LL1 = (TempIn * (1 - Tc / Tb) - x3LL1)/Tb;
    TempIn = TempIn * Tc / Tb + x3LL1; // output of first lead lag
  }
  // State 4
  if (Tb1 < TS_THRESHOLD * t_inc) {
//    printf ("entering Tb < 4h\n");
    dx4LL2 = 0;
    TempIn = TempIn; // output of second lead lag - just pass input
  } else {
    dx4LL2 = (TempIn * (1 - Tc1 / Tb1) - x4LL2)/Tb1;
    TempIn = TempIn * Tc1 / Tb1 + x4LL2; // output of second lead lag
  }
  // State 1
  if (Ta < TS_THRESHOLD * t_inc) x1Va = Ka * TempIn;
  if (x1Va > Vamax) x1Va = Vamax; 
  if (x1Va < Vamin) x1Va = Vamin; 
  if (Ta < TS_THRESHOLD * t_inc) dx1Va = 0;
  else {
     dx1Va = (Ka * TempIn - x1Va) / Ta;
  }
  if (dx1Va > 0 && x1Va >= Vamax) dx1Va = 0;
  if (dx1Va < 0 && x1Va <= Vamin) dx1Va = 0;

  x1Va_1 = x1Va + dx1Va * t_inc;
  x2Vcomp_1 = x2Vcomp + dx2Vcomp * t_inc;
  x3LL1_1 = x3LL1 + dx3LL1 * t_inc;
  x4LL2_1 = x4LL2 + dx4LL2 * t_inc;
  x5Deriv_1 = x5Deriv + dx5Deriv * t_inc;

  printf("esst1a dx: %f\t%f\t%f\t%f\t%f\n", dx1Va, dx2Vcomp, dx3LL1, dx4LL2, dx5Deriv);
  printf("esst1a x: %f\t%f\t%f\t%f\t%f\n", x1Va_1, x2Vcomp_1, x3LL1_1, x4LL2_1, x5Deriv_1);

  double Temp = Klr * (LadIfd - Ilr);
  if (Temp < 0) Temp = 0;
  Temp = x1Va - Temp;
  // Ignore Over and Under Excitation Limit for now
  if (Temp > (Vterm * Vrmax - Kc * LadIfd)) Temp = Vterm * Vrmax - Kc * LadIfd;
  if (Temp < (Vterm * Vrmin)) Temp = Vterm * Vrmin;
  Efd = Temp;

  //printf("esst1a Efd: %f\n", Efd);*/
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Esst1aModel::corrector(double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR);
 /* // State 2
  if (Tr < TS_THRESHOLD * t_inc) {
    x2Vcomp_1 = Vcomp; // Must propogate input value instantaneously
    dx2Vcomp_1 = 0;
  } else {
    dx2Vcomp_1 = 1 / Tr * (Vcomp - x2Vcomp_1);
  }
  // State 5
  double TempIn;
  if (Kf > 0) {
    //TempIn = Klr * (LadIfd - Ilr);
    //if (TempIn < 0) TempIn = 0;
    // Note: if Ta = 0, then we would have an algebriac loop here which
    // could cause numerical troubles. Thus we enforce Ta > 0 if (Kf > 0 and Tf > 0)
    //TempIn = x1Va_1 - TempIn;
    //if ("VOS at Output") TempIn = TempIn + Vstab; // TBD: "VOS at Output"???
    // Ignore Over and Under Excitation Limit for now
	
	TempIn = x1Va_1;
    double UseTf;
    if (Tf > TS_THRESHOLD * t_inc) UseTf = Tf;
    else UseTf = TS_THRESHOLD * t_inc;
    dx5Deriv_1 = (TempIn * (-Kf / UseTf) - x5Deriv_1) / UseTf;
    TempIn = TempIn * Kf / UseTf + x5Deriv_1; // Input into summation block
  } else {
    dx5Deriv_1 = 0;
    TempIn = 0; // Input into summation block;
  }
  // State 3
  TempIn = - x2Vcomp_1 - TempIn + Vref;
  TempIn = TempIn + Vstab;
  //if ("VOS at Input") TempIn = TempIn + Vstab; // TBD: "VOS at Input"???
  if (TempIn > Vimax) TempIn = Vimax;
  if (TempIn < Vimin) TempIn = Vimin;
  // Ignore Under Excitation Limit
  if (Tb < TS_THRESHOLD * t_inc) {
    dx3LL1_1 = 0;
    TempIn = TempIn; // output of first lead lag - just pass input
  } else {
    dx3LL1_1 = (TempIn * (1 - Tc / Tb) - x3LL1_1)/Tb;
    TempIn = TempIn * Tc / Tb + x3LL1_1; // output of first lead lag
  }
  // State 4
  if (Tb1 < TS_THRESHOLD * t_inc) {
    dx4LL2_1 = 0;
    TempIn = TempIn; // output of second lead lag - just pass input
  } else {
    dx4LL2_1 = (TempIn * (1 - Tc1 / Tb1) - x4LL2_1)/Tb1;
    TempIn = TempIn * Tc1 / Tb1 + x4LL2_1; // output of second lead lag
  }
  // State 1
  if (Ta < TS_THRESHOLD * t_inc) x1Va_1 = Ka * TempIn;
  if (x1Va_1 > Vamax) x1Va_1 = Vamax; 
  if (x1Va_1 < Vamin) x1Va_1 = Vamin; 
  if (Ta < TS_THRESHOLD * t_inc) dx1Va_1 = 0;
  else dx1Va_1 = (Ka * TempIn - x1Va_1) / Ta;
  if (dx1Va_1 > 0 && x1Va_1 >= Vamax) dx1Va_1 = 0;
  if (dx1Va_1 < 0 && x1Va_1 <= Vamin) dx1Va_1 = 0;


  x1Va_1 = x1Va + (dx1Va + dx1Va_1) / 2.0 * t_inc;
  x2Vcomp_1 = x2Vcomp + (dx2Vcomp + dx2Vcomp_1) / 2.0 * t_inc;
  x3LL1_1 = x3LL1 + (dx3LL1 + dx3LL1_1) / 2.0 * t_inc;
  x4LL2_1 = x4LL2 + (dx4LL2 + dx4LL2_1) / 2.0 * t_inc;
  x5Deriv_1 = x5Deriv + (dx5Deriv + dx5Deriv_1) / 2.0 * t_inc;

  //printf("esst1a dx: %f\t%f\t%f\t%f\t%f\n", dx1Va_1, dx2Vcomp_1, dx3LL1_1, dx4LL2_1, dx5Deriv_1);
//  printf("esst1a x: %f\t%f\t%f\t%f\t%f\n", x1Va_1, x2Vcomp_1, x3LL1_1, x4LL2_1, x5Deriv_1);
  
  double Temp = Klr * (LadIfd - Ilr);
  if (Temp < 0) Temp = 0;
  Temp = x1Va_1 - Temp;
  // Ignore Over and Under Excitation Limit for now
  if (Temp > (Vterm * Vrmax - Kc * LadIfd)) Temp = Vterm * Vrmax - Kc * LadIfd;
  if (Temp < (Vterm * Vrmin)) Temp = Vterm * Vrmin;
  Efd = Temp;

  //printf("esst1a Efd: %f\n", Efd);*/
}

/**
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void gridpack::dynamic_simulation::Esst1aModel::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

/**
 * Set the field current parameter inside the exciter
 * @param fldc value of the field current
 */
void gridpack::dynamic_simulation::Esst1aModel::setFieldCurrent(double fldc)
{
  LadIfd = fldc;
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double gridpack::dynamic_simulation::Esst1aModel::getFieldVoltage()
{
  return Efd;
}

/** 
 * Get the value of the field current parameter
 * @return value of field current
 */
double gridpack::dynamic_simulation::Esst1aModel::getFieldCurrent()
{
  return 0.0;
}

/** 
 * Set the value of the Vterminal
 * @return value of field current
 */
void gridpack::dynamic_simulation::Esst1aModel::setVterminal(double mag)
{
  //Vterminal = mag;
  Vterm = mag;
}

/** 
 * Set the value of the omega
 * @return value of field current
 */
void gridpack::dynamic_simulation::Esst1aModel::setOmega(double omega)
{
  //w = omega;
}

/** 
 * Set the value of the Vcomp
 * @return value of the Vcomp
 */
void gridpack::dynamic_simulation::Esst1aModel::setVcomp(double vtmp)
{
  Vcomp = vtmp;
}

void gridpack::dynamic_simulation::Esst1aModel::setVstab(double vtmp)
{
  Vstab = vtmp;
}

void gridpack::dynamic_simulation::Esst1aModel::setVothsg(double vtmp)
{
  Vothsg = vtmp;
}

void gridpack::dynamic_simulation::Esst1aModel::setVuel(double vtmp)
{
  Vuel = vtmp;
}

void gridpack::dynamic_simulation::Esst1aModel::setVoel(double vtmp)
{
  Voel = vtmp;
}

/**
 * Set internal state parameter in exciter
 * @param name character string corresponding to state variable
 * @param value new value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::Esst1aModel::setState(std::string name,
    double value)
{
  return false;
}

/**
 * Get internal state parameter in exciter
 * @param name character string corresponding to state variable
 * @param value current value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::Esst1aModel::getState(std::string name,
    double *value)
{
  return false;
}
