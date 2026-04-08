/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   esdc2a.cpp
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief  ESDC2A exciter model (IEEE DC2A / PSS/E ESDC2A).
 *
 * Identical topology to EXDC2.  The only difference is that the
 * voltage-regulator anti-windup limits scale with terminal voltage Vt:
 *   VRmax_eff = Vrmax_c * Vt   (Vrmax_c = VRMAX if VRMAX != 0, else 999)
 *   VRmin_eff = VRMIN * Vt
 *
 * This matches the ANDES ESDC2A implementation (esdc2a.py):
 *   VRU = VarService('VRMAXc * v')
 *   VRL = VarService('VRMIN * v')
 *
 * PSS/E autocorrections applied:
 *   - Swap Vrmax/Vrmin if inverted
 *   - VRMAX == 0 treated as unlimited (999 pu)
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <algorithm>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_exciter_model.hpp"
#include "esdc2a.hpp"

#define TS_THRESHOLD 4

gridpack::dynamic_simulation::Esdc2aModel::Esdc2aModel(void)
{
  Vs = 0.0;
  sat_A = 0.0;
  sat_B = 0.0;
  has_Sat = true;
  has_leadlag = true;
  zero_TA = false;
  zero_TR = false;
}

gridpack::dynamic_simulation::Esdc2aModel::~Esdc2aModel(void)
{
}

void gridpack::dynamic_simulation::Esdc2aModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(EXCITER_TR,    &TR,          idx)) TR          = 0.0;
  if (!data->getValue(EXCITER_KA,    &KA,          idx)) KA          = 0.0;
  if (!data->getValue(EXCITER_TA,    &TA,          idx)) TA          = 0.0;
  if (!data->getValue(EXCITER_TB,    &TB,          idx)) TB          = 0.0;
  if (!data->getValue(EXCITER_TC,    &TC,          idx)) TC          = 0.0;
  if (!data->getValue(EXCITER_VRMAX, &Vrmax_param, idx)) Vrmax_param = 0.0;
  if (!data->getValue(EXCITER_VRMIN, &Vrmin_param, idx)) Vrmin_param = 0.0;
  if (!data->getValue(EXCITER_KE,    &KE,          idx)) KE          = 0.0;
  if (!data->getValue(EXCITER_TE,    &TE,          idx)) TE          = 0.0;
  if (!data->getValue(EXCITER_KF,    &KF,          idx)) KF          = 0.0;
  if (!data->getValue(EXCITER_TF1,   &TF1,         idx)) TF1         = 0.0;
  if (!data->getValue(EXCITER_E1,    &E1,          idx)) E1          = 0.0;
  if (!data->getValue(EXCITER_SE1,   &SE1,         idx)) SE1         = 0.0;
  if (!data->getValue(EXCITER_E2,    &E2,          idx)) E2          = 0.0;
  if (!data->getValue(EXCITER_SE2,   &SE2,         idx)) SE2         = 0.0;

  // VRMAX == 0 means unlimited in PSS/E ESDC2A convention
  Vrmax_c = (fabs(Vrmax_param) < 1e-6) ? 999.0 : Vrmax_param;

  // Autocorrection: swap if inverted
  if (Vrmax_c < Vrmin_param) {
    double tmp = Vrmax_c; Vrmax_c = Vrmin_param; Vrmin_param = tmp;
  }

  if (fabs(SE1 * SE2) < 1e-6) has_Sat = false;
  if (fabs(TB)        < 1e-6) has_leadlag = false;
  if (fabs(TR)        < 1e-6) zero_TR = true;
  if (fabs(TA)        < 1e-6) zero_TA = true;

  if (has_Sat) {
    double r = std::sqrt(SE1 * E1 / (SE2 * E2));
    sat_A = (E1 - r * E2) / (1.0 - r);
    sat_B = SE1 * E1 / ((E1 - sat_A) * (E1 - sat_A));
  }

  // Set up blocks with nominal limits (will be Vt-scaled at runtime)
  if (!zero_TR) {
    Vmeas_blk.setparams(1.0, TR);
  }
  if (has_leadlag) {
    Leadlag_blk.setparams(TC, TB);
  }
  // Regulator limits at nominal Vt=1 pu; updated each step in computeModel
  if (!zero_TA) {
    Regulator_blk.setparams(KA, TA, Vrmin_param, Vrmax_c, -1000.0, 1000.0);
  } else {
    Regulator_gain_blk.setparams(KA, Vrmin_param, Vrmax_c);
  }

  double a[2], b[2];
  a[0] = TF1; a[1] = 1.0;
  b[0] = KF;  b[1] = 0.0;
  Feedback_blk.setcoeffs(a, b);

  Output_blk.setparams(TE);
}

double gridpack::dynamic_simulation::Esdc2aModel::Sat(double x)
{
  if (x < sat_A) return 0.0;
  return sat_B * (x - sat_A) * (x - sat_A) / x;
}

void gridpack::dynamic_simulation::Esdc2aModel::init(double Vm, double Va, double ts)
{
  double mult_ts = TS_THRESHOLD * ts;

  if (TR > 0.0 && TR < 0.25 * mult_ts) {
    TR = 0.0; zero_TR = true;
  } else if (TR > 0.25 * mult_ts && TR < 0.5 * mult_ts) {
    TR = 0.5 * mult_ts;
    Vmeas_blk.setparams(1.0, TR);
  } else if (fabs(TR) < 1e-6) {
    zero_TR = true;
  }

  if (TB > 0.0 && TB < 0.5 * mult_ts) {
    TB = 0.0; has_leadlag = false;
  } else if (TB > 0.5 * mult_ts && TB < mult_ts) {
    TB = mult_ts;
    if (has_leadlag) Leadlag_blk.setparams(TC, TB);
  }

  if (TE > 0.0 && TE < mult_ts) {
    TE = mult_ts;
    Output_blk.setparams(TE);
  }

  if (TF1 > 0.0 && TF1 < mult_ts) {
    TF1 = mult_ts;
    double a[2], b[2];
    a[0] = TF1; a[1] = 1.0;
    b[0] = KF;  b[1] = 0.0;
    Feedback_blk.setcoeffs(a, b);
  }

  if (TA > 0.0 && TA < mult_ts) {
    TA = mult_ts; zero_TA = false;
  } else if (fabs(TA) < 1e-6) {
    zero_TA = true;
  }

  // Vt-scaled limits at initialization (Ec already set to Vt by setVterminal)
  double vrmax_eff = Vrmax_c  * Ec;
  double vrmin_eff = Vrmin_param * Ec;
  if (!zero_TA) {
    Regulator_blk.setparams(KA, TA, vrmin_eff, vrmax_eff, -1000.0, 1000.0);
  } else {
    Regulator_gain_blk.setparams(KA, vrmin_eff, vrmax_eff);
  }

  VF = Feedback_blk.init_given_u(Efd);

  double output_blk_in = Output_blk.init_given_y(Efd);

  double sat_signal = KE * Efd;
  if (has_Sat) sat_signal += Efd * Sat(Efd);

  VR = output_blk_in + sat_signal;

  // Widen effective limits if initial VR is outside bounds
  if (VR > vrmax_eff) {
    vrmax_eff = VR;
    if (!zero_TA) Regulator_blk.setparams(KA, TA, vrmin_eff, vrmax_eff, -1000.0, 1000.0);
    else Regulator_gain_blk.setparams(KA, vrmin_eff, vrmax_eff);
  }
  if (VR < vrmin_eff) {
    vrmin_eff = VR;
    if (!zero_TA) Regulator_blk.setparams(KA, TA, vrmin_eff, vrmax_eff, -1000.0, 1000.0);
    else Regulator_gain_blk.setparams(KA, vrmin_eff, vrmax_eff);
  }

  if (zero_TA) {
    VLL = VR / KA;
  } else {
    VLL = Regulator_blk.init_given_y(VR);
  }

  double leadlag_blk_in;
  if (has_leadlag) {
    leadlag_blk_in = Leadlag_blk.init_given_y(VLL);
  } else {
    leadlag_blk_in = VLL;
  }

  double Verr = leadlag_blk_in - VF + Vs;

  if (!zero_TR) {
    Vmeas = Vmeas_blk.init_given_u(Ec);
  } else {
    Vmeas = Ec;
  }

  Vref = Verr + Vmeas;
}

void gridpack::dynamic_simulation::Esdc2aModel::computeModel(
    double t_inc, IntegrationStage int_flag)
{
  if (!zero_TR) {
    Vmeas = Vmeas_blk.getoutput(Ec, t_inc, int_flag, true);
  } else {
    Vmeas = Ec;
  }

  double Verr = Vref - Vmeas + Vs;

  VF = Feedback_blk.getoutput(Efd, t_inc, int_flag, true);

  double leadlag_blk_in = Verr - VF;
  if (has_leadlag) {
    VLL = Leadlag_blk.getoutput(leadlag_blk_in, t_inc, int_flag, true);
  } else {
    VLL = leadlag_blk_in;
  }

  // Update Vt-scaled regulator limits each step
  double vrmax_eff = Vrmax_c    * Ec;
  double vrmin_eff = Vrmin_param * Ec;
  if (zero_TA) {
    Regulator_gain_blk.setparams(KA, vrmin_eff, vrmax_eff);
    VR = Regulator_gain_blk.getoutput(VLL);
  } else {
    Regulator_blk.setxlimits(vrmin_eff, vrmax_eff);
    VR = Regulator_blk.getoutput(VLL, t_inc, int_flag, true);
  }

  double sat_signal = KE * Efd;
  if (has_Sat) sat_signal += Efd * Sat(Efd);

  double output_blk_in = VR - sat_signal;

  Efd = Output_blk.getoutput(output_blk_in, t_inc, int_flag, true);
}

void gridpack::dynamic_simulation::Esdc2aModel::predictor(double t_inc, bool flag)
{
  computeModel(t_inc, PREDICTOR);
}

void gridpack::dynamic_simulation::Esdc2aModel::corrector(double t_inc, bool flag)
{
  computeModel(t_inc, CORRECTOR);
}

void gridpack::dynamic_simulation::Esdc2aModel::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

double gridpack::dynamic_simulation::Esdc2aModel::getFieldVoltage()
{
  return Efd;
}

void gridpack::dynamic_simulation::Esdc2aModel::setVterminal(double Vm)
{
  Ec = Vm;
}

void gridpack::dynamic_simulation::Esdc2aModel::setVstab(double Vstab)
{
  Vs = Vstab;
}

void gridpack::dynamic_simulation::Esdc2aModel::setExtBusNum(int ExtBusNum)
{
  p_bus_num = ExtBusNum;
}

void gridpack::dynamic_simulation::Esdc2aModel::setExtGenId(std::string ExtGenId)
{
  p_gen_id = ExtGenId;
}

bool gridpack::dynamic_simulation::Esdc2aModel::setState(std::string name,
    double value)
{
  return false;
}

bool gridpack::dynamic_simulation::Esdc2aModel::getState(std::string name,
    double *value)
{
  return false;
}
