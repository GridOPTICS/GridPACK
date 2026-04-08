/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   exst1.cpp
 * @author Yousu Chen
 * @Updated:  March 2026
 *
 * @brief  EXST1 static exciter model.
 *
 * Signal flow:
 *   Vmeas  = 1/(1+sTR) * Vt
 *   Verr   = Vref + Vs - Vmeas - VF
 *   VI     = clamp(Verr, VIMIN, VIMAX)
 *   VLL    = (1+sTC)/(1+sTB) * VI      [bypassed if TB = 0]
 *   VR     = KA/(1+sTA) * VLL          [anti-windup VRMIN, VRMAX]
 *   Efd    = VR - KC * Ifd
 *   VF     = KF * s/(1+sTF) * Efd      [rate feedback]
 *
 * PSS/E autocorrections:
 *   - Swap VRMAX/VRMIN if inverted
 *   - Swap VIMAX/VIMIN if inverted
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
#include "exst1.hpp"

#define TS_THRESHOLD 4

gridpack::dynamic_simulation::Exst1Model::Exst1Model(void)
{
  Vs = 0.0;
  LadIfd = 0.0;
  zero_TR = false;
  zero_TB = false;
  zero_TA = false;
}

gridpack::dynamic_simulation::Exst1Model::~Exst1Model(void)
{
}

void gridpack::dynamic_simulation::Exst1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(EXCITER_TR,    &TR,    idx)) TR    = 0.0;
  if (!data->getValue(EXCITER_VIMAX, &VIMAX, idx)) VIMAX = 1.0;
  if (!data->getValue(EXCITER_VIMIN, &VIMIN, idx)) VIMIN = -1.0;
  if (!data->getValue(EXCITER_TC,    &TC,    idx)) TC    = 0.0;
  if (!data->getValue(EXCITER_TB,    &TB,    idx)) TB    = 0.0;
  if (!data->getValue(EXCITER_KA,    &KA,    idx)) KA    = 0.0;
  if (!data->getValue(EXCITER_TA,    &TA,    idx)) TA    = 0.0;
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 0.0;
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = 0.0;
  if (!data->getValue(EXCITER_KC,    &KC,    idx)) KC    = 0.0;
  if (!data->getValue(EXCITER_KF,    &KF,    idx)) KF    = 0.0;
  if (!data->getValue(EXCITER_TF1,   &TF,    idx)) TF    = 0.0;

  // Autocorrection: swap limits if inverted
  if (Vrmax < Vrmin) { double tmp = Vrmax; Vrmax = Vrmin; Vrmin = tmp; }
  if (VIMAX < VIMIN) { double tmp = VIMAX; VIMAX = VIMIN; VIMIN = tmp; }

  if (fabs(TR) < 1e-6) zero_TR = true;
  if (fabs(TB) < 1e-6) zero_TB = true;
  if (fabs(TA) < 1e-6) zero_TA = true;

  // Set up blocks
  if (!zero_TR) {
    Vmeas_blk.setparams(1.0, TR);
  }
  if (!zero_TB) {
    Leadlag_blk.setparams(TC, TB);
  }
  if (!zero_TA) {
    Regulator_blk.setparams(KA, TA, Vrmin, Vrmax, -1000.0, 1000.0);
  } else {
    Regulator_gain_blk.setparams(KA, Vrmin, Vrmax);
  }
  double a[2], b[2];
  a[0] = TF; a[1] = 1.0;
  b[0] = KF; b[1] = 0.0;
  Feedback_blk.setcoeffs(a, b);
}

void gridpack::dynamic_simulation::Exst1Model::init(double Vm, double Va, double ts)
{
  double mult_ts = TS_THRESHOLD * ts;

  // TR
  if (TR > 0.0 && TR < 0.25 * mult_ts) {
    TR = 0.0; zero_TR = true;
  } else if (TR > 0.25 * mult_ts && TR < 0.5 * mult_ts) {
    TR = 0.5 * mult_ts;
    Vmeas_blk.setparams(1.0, TR);
  } else if (fabs(TR) < 1e-6) {
    zero_TR = true;
  }

  // TB
  if (TB > 0.0 && TB < 0.5 * mult_ts) {
    TB = 0.0; zero_TB = true;
  } else if (TB > 0.5 * mult_ts && TB < mult_ts) {
    TB = mult_ts;
    if (!zero_TB) Leadlag_blk.setparams(TC, TB);
  }

  // TF
  if (TF > 0.0 && TF < mult_ts) {
    TF = mult_ts;
    double a[2], b[2];
    a[0] = TF; a[1] = 1.0;
    b[0] = KF; b[1] = 0.0;
    Feedback_blk.setcoeffs(a, b);
  }

  // TA
  if (TA > 0.0 && TA < mult_ts) {
    TA = mult_ts; zero_TA = false;
    Regulator_blk.setparams(KA, TA, Vrmin, Vrmax, -1000.0, 1000.0);
  } else if (fabs(TA) < 1e-6) {
    zero_TA = true;
  }

  // Initialize feedback block from Efd
  VF = Feedback_blk.init_given_u(Efd);

  // VR = Efd + KC * Ifd  (no output integrator)
  VR = Efd + KC * LadIfd;

  // Widen VR limits if needed
  if (VR > Vrmax) {
    Vrmax = VR;
    if (!zero_TA) Regulator_blk.setparams(KA, TA, Vrmin, Vrmax, -1000.0, 1000.0);
    else Regulator_gain_blk.setparams(KA, Vrmin, Vrmax);
  }
  if (VR < Vrmin) {
    Vrmin = VR;
    if (!zero_TA) Regulator_blk.setparams(KA, TA, Vrmin, Vrmax, -1000.0, 1000.0);
    else Regulator_gain_blk.setparams(KA, Vrmin, Vrmax);
  }

  // Initialize regulator from VR
  double leadlag_in;
  if (zero_TA) {
    leadlag_in = VR / KA;
  } else {
    leadlag_in = Regulator_blk.init_given_y(VR);
  }

  // Initialize lead-lag from its output (=leadlag_in since steady-state)
  double vi;
  if (!zero_TB) {
    vi = Leadlag_blk.init_given_y(leadlag_in);
  } else {
    vi = leadlag_in;
  }

  // VI = clamp(Verr, VIMIN, VIMAX)  → Verr = vi at steady state
  double Verr = vi;
  // Widen VI limits if initial error violates them
  if (Verr > VIMAX) VIMAX = Verr;
  if (Verr < VIMIN) VIMIN = Verr;

  // Vmeas
  if (!zero_TR) {
    Vmeas = Vmeas_blk.init_given_u(Ec);
  } else {
    Vmeas = Ec;
  }

  Vref = Verr + Vmeas + VF - Vs;
}

void gridpack::dynamic_simulation::Exst1Model::computeModel(
    double t_inc, IntegrationStage int_flag)
{
  // Voltage measurement
  if (!zero_TR) {
    Vmeas = Vmeas_blk.getoutput(Ec, t_inc, int_flag, true);
  } else {
    Vmeas = Ec;
  }

  // Rate feedback
  VF = Feedback_blk.getoutput(Efd, t_inc, int_flag, true);

  // Error signal + input limiter
  double Verr = Vref + Vs - Vmeas - VF;
  double VI = std::max(VIMIN, std::min(VIMAX, Verr));

  // Lead-lag
  if (!zero_TB) {
    VLL = Leadlag_blk.getoutput(VI, t_inc, int_flag, true);
  } else {
    VLL = VI;
  }

  // Regulator
  if (zero_TA) {
    VR = Regulator_gain_blk.getoutput(VLL);
  } else {
    VR = Regulator_blk.getoutput(VLL, t_inc, int_flag, true);
  }

  // Output: no integrator, subtract KC * Ifd
  Efd = VR - KC * LadIfd;
}

void gridpack::dynamic_simulation::Exst1Model::predictor(double t_inc, bool flag)
{
  computeModel(t_inc, PREDICTOR);
}

void gridpack::dynamic_simulation::Exst1Model::corrector(double t_inc, bool flag)
{
  computeModel(t_inc, CORRECTOR);
}

void gridpack::dynamic_simulation::Exst1Model::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

double gridpack::dynamic_simulation::Exst1Model::getFieldVoltage()
{
  return Efd;
}

void gridpack::dynamic_simulation::Exst1Model::setFieldCurrent(double fldc)
{
  LadIfd = fldc;
}

void gridpack::dynamic_simulation::Exst1Model::setVterminal(double Vm)
{
  Ec = Vm;
}

void gridpack::dynamic_simulation::Exst1Model::setVstab(double Vstab)
{
  Vs = Vstab;
}

void gridpack::dynamic_simulation::Exst1Model::setExtBusNum(int ExtBusNum)
{
  p_bus_num = ExtBusNum;
}

void gridpack::dynamic_simulation::Exst1Model::setExtGenId(std::string ExtGenId)
{
  p_gen_id = ExtGenId;
}

bool gridpack::dynamic_simulation::Exst1Model::setState(std::string name,
    double value)
{
  return false;
}

bool gridpack::dynamic_simulation::Exst1Model::getState(std::string name,
    double *value)
{
  return false;
}
