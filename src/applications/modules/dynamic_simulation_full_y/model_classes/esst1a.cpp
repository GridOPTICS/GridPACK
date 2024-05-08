/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   esst1a.cpp
 *
 * @brief  ESST1A model
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

#define TS_THRESHOLD 4

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Esst1aModel::Esst1aModel(void)
{
  zero_TF = false;
  zero_TB = false;
  zero_TB1 = false;
  OptionToModifyLimitsForInitialStateLimitViolation = true;

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

  // right now we just hard code UEL, VOS, Vuel, Voel and Vothsg(Vstab)
  Vothsg = 0.0;
  UEL = 1.0;
  VOS = 1.0;
  Vuel = 0.0;
  Voel = 1000.0;

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
	
	return 0.0;
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
  if (Tf < TS_THRESHOLD * ts) zero_TF = true;
  if (Tb < TS_THRESHOLD * ts) zero_TB = true;
  if (Tb1 < TS_THRESHOLD * ts) zero_TB1 = true;
  if (Ta < TS_THRESHOLD * ts) zero_TA = true;
  if (Tr < TS_THRESHOLD * ts) zero_TR = true;
  
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
  
  Vterm = mag;
  
  Vf = Feedback_blkF.init_given_u(Efd);
  
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (Efd > (Vterm * Vrmax - Kc * LadIfd)) Vrmax = (Efd + Kc * LadIfd) / Vterm+0.21;
    if (Efd < (Vterm * Vrmin)) Vrmin = Efd / Vterm-0.1;
  }
  
  // LV Gate?
  // assume LV gate always take feedforward path during initialization, may need to adjust Voel
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
	  if (Efd > Voel) Voel = Efd + 0.1;
	  LVGate_blk.setparams(Voel);
  }

  // HV Gate?
  // assume HV gate always take feedforward path during initialization, may need to adjust Vuel
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
	  if (Efd < Vuel) Vuel = Efd - 0.1;
	  HVGate_blk2.setparams(Vuel);
  }

  if (VOS==2.0) {
	  if ((LadIfd - Ilr) * Klr > 0.0) VA = (LadIfd - Ilr) * Klr - Vothsg + Efd;
	  else VA = - Vothsg + Efd;
  } else {
	  if ((LadIfd - Ilr) * Klr > 0.0) VA = (LadIfd - Ilr) * Klr + Efd;
	  else VA = Efd;
  }


  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (VA > Vamax) Vamax = VA+0.1;
    if (VA < Vamin) Vamin = VA-0.1;
  }
  VLL1 = Regulator_blk.init_given_y(VA);


  VLL = Leadlag_blkBC1.init_given_y(VLL1);
  double u1 = Leadlag_blkBC.init_given_y(VLL);

  // HV Gate?
  // assume HV gate always take feedforward path during initialization, may need to adjust Vuel
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
	  if (u1 < Vuel) Vuel = u1 - 0.1;
	  HVGate_blk1.setparams(Vuel);
  }
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (u1 > Vimax) Vimax = u1+0.1;
    if (u1 < Vimin) Vimin = u1-0.1;
  }

  double Vop = 0.0;
  if (UEL == 1.0) Vop += Vuel;
  if (VOS == 1.0) Vop += Vothsg;
  Vmeas = Filter_blkR.init_given_u(Vcomp);
  Vref = u1 + Vmeas - Vop + Vf;

}

/**
 * computeModel - Updates the model states and gets the output
 */
void gridpack::dynamic_simulation::Esst1aModel::computeModel(double t_inc,IntegrationStage int_flag)
{
  
  if(!zero_TR) {
    Vmeas = Filter_blkR.getoutput(Vcomp, t_inc, int_flag, true);
  } else {
    Vmeas = Vcomp;
  }
  double Vop = 0.0;
  if (UEL == 1.0) Vop += Vuel;
  if (VOS == 1.0) Vop += Vothsg;
  double Verr = Vref - Vmeas + Vop;

  Vf = Feedback_blkF.getoutput(Efd, t_inc, int_flag, true);

  double leadlag_blk_in = Verr - Vf;

  if (leadlag_blk_in > Vimax)
    leadlag_blk_in = Vimax;
  else if (leadlag_blk_in < Vimin)
    leadlag_blk_in = Vimin;

  if (UEL == 2.0) leadlag_blk_in = HVGate_blk1.getoutput(leadlag_blk_in);

  VLL = Leadlag_blkBC.getoutput(leadlag_blk_in, t_inc, int_flag, true);
  
  VLL1 = Leadlag_blkBC1.getoutput(VLL, t_inc, int_flag, true); 

  if (zero_TA) {
    VA = Regulator_gain_blk.getoutput(VLL1);
  } else {
    VA = Regulator_blk.getoutput(VLL, t_inc, int_flag, true);
  }

  double u1 = 0.0;
  if (VOS==2.0) {
    if ((LadIfd - Ilr) * Klr > 0.0) u1 = VA + Vothsg - (LadIfd - Ilr) * Klr; 
    else u1 = VA + Vothsg; 
  } else {
    if ((LadIfd - Ilr) * Klr > 0.0) u1 = VA - (LadIfd - Ilr) * Klr; 
    else u1 = VA; 
  }

  double u2 = 0.0;
  if (UEL == 3.0) u2 = HVGate_blk2.getoutput(u1);
  else u2 = u1;

  u2 = LVGate_blk.getoutput(u2);

  double VT = Vterm;
  
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

}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Esst1aModel::corrector(double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR);

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
