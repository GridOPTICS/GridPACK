/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   ieeet2.cpp
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief  IEEET2 — IEEE Type 2 DC Exciter
 *
 *  Same topology as IEEET1 with two differences:
 *   1. Rate feedback has an additional first-order lag 1/(1+TF2*s) in cascade:
 *        Efd → 1/(1+TF2*s) → KF*TF1*s/(1+TF1*s) → VF
 *   2. VR output limits scale with terminal voltage Vt:
 *        VRmax_eff = VRMAX * Vt,  VRmin_eff = VRMIN * Vt
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <cmath>
#include <algorithm>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_exciter_model.hpp"
#include "ieeet2.hpp"

#define TS_THRESHOLD 4

gridpack::dynamic_simulation::Ieeet2Model::Ieeet2Model(void)
{
  Vs = 0.0;
  sat_A = 0.0;
  sat_B = 0.0;
  has_Sat = true;
  zero_TA  = false;
  zero_TR  = false;
  zero_TF2 = false;
}

gridpack::dynamic_simulation::Ieeet2Model::~Ieeet2Model(void)
{
}

void gridpack::dynamic_simulation::Ieeet2Model::load(
    boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  if (!data->getValue(EXCITER_TR,    &TR,    idx)) TR    = 0.0;
  if (!data->getValue(EXCITER_KA,    &KA,    idx)) KA    = 0.0;
  if (!data->getValue(EXCITER_TA,    &TA,    idx)) TA    = 0.0;
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 0.0;
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = 0.0;
  if (!data->getValue(EXCITER_KE,    &KE,    idx)) KE    = 0.0;
  if (!data->getValue(EXCITER_TE,    &TE,    idx)) TE    = 0.0;
  if (!data->getValue(EXCITER_KF,    &KF,    idx)) KF    = 0.0;
  if (!data->getValue(EXCITER_TF1,   &TF1,   idx)) TF1   = 0.0;
  if (!data->getValue(EXCITER_TF2,   &TF2,   idx)) TF2   = 0.0;
  if (!data->getValue(EXCITER_E1,    &E1,    idx)) E1    = 0.0;
  if (!data->getValue(EXCITER_SE1,   &SE1,   idx)) SE1   = 0.0;
  if (!data->getValue(EXCITER_E2,    &E2,    idx)) E2    = 0.0;
  if (!data->getValue(EXCITER_SE2,   &SE2,   idx)) SE2   = 0.0;

  if (fabs(SE1*SE2) < 1e-6) has_Sat = false;
  if (fabs(TA)  < 1e-6) zero_TA  = true;
  if (fabs(TR)  < 1e-6) zero_TR  = true;
  if (fabs(TF2) < 1e-6) zero_TF2 = true;

  if (has_Sat) {
    double r = std::sqrt(SE1*E1 / (SE2*E2));
    sat_A = (E1 - r*E2) / (1.0 - r);
    sat_B = SE1*E1 / ((E1 - sat_A)*(E1 - sat_A));
  }

  if (!zero_TR)  Vmeas_blk.setparams(1.0, TR);
  if (!zero_TA)  Regulator_blk.setparams(KA, TA, Vrmin, Vrmax, -1000.0, 1000.0);
  else           Regulator_gain_blk.setparams(KA, Vrmin, Vrmax);
  if (!zero_TF2) TF2_blk.setparams(1.0, TF2);

  double a[2], b[2];
  a[0] = TF1; a[1] = 1.0;
  b[0] = KF;  b[1] = 0.0;
  Feedback_blk.setcoeffs(a, b);

  Output_blk.setparams(TE);
}

double gridpack::dynamic_simulation::Ieeet2Model::Sat(double x)
{
  if (x < sat_A) return 0.0;
  return sat_B*(x - sat_A)*(x - sat_A)/x;
}

void gridpack::dynamic_simulation::Ieeet2Model::init(
    double Vm, double Va, double ts)
{
  // Inner feedback lag: at steady state TF2_lag output = Efd
  // Outer washout: at steady state VF = 0
  VF = Feedback_blk.init_given_u(Efd);
  if (!zero_TF2) TF2_blk.init_given_u(Efd);

  double output_blk_in = Output_blk.init_given_y(Efd);
  double sat_signal = KE*Efd;
  if (has_Sat) sat_signal += Efd*Sat(Efd);

  VR = output_blk_in + sat_signal;

  // At init Vt≈1.0, so effective limits ≈ fixed limits
  double vrmax_eff = Vrmax * Vm;
  double vrmin_eff = Vrmin * Vm;
  VR = std::min(vrmax_eff, std::max(vrmin_eff, VR));

  double Regulator_blk_in;
  if (zero_TA) {
    Regulator_blk_in = VR / KA;
  } else {
    Regulator_blk_in = Regulator_blk.init_given_y(VR);
  }

  double Verr = Regulator_blk_in + VF;

  if (!zero_TR) {
    Vmeas = Vmeas_blk.init_given_u(Ec);
  } else {
    Vmeas = Ec;
  }

  Vref = Verr + Vmeas - Vs;
}

void gridpack::dynamic_simulation::Ieeet2Model::computeModel(
    double t_inc, IntegrationStage int_flag)
{
  if (!zero_TR) {
    Vmeas = Vmeas_blk.getoutput(Ec, t_inc, int_flag, true);
  } else {
    Vmeas = Ec;
  }

  double Verr = Vref - Vmeas + Vs;

  // Two-stage rate feedback: Efd → TF2 lag → TF1 washout → VF
  double vf_inner;
  if (!zero_TF2) {
    vf_inner = TF2_blk.getoutput(Efd, t_inc, int_flag, true);
  } else {
    vf_inner = Efd;
  }
  VF = Feedback_blk.getoutput(vf_inner, t_inc, int_flag, true);

  double Regulator_blk_in = Verr - VF;

  // VR limits scale with Vt (key IEEET2 difference)
  double vrmax_eff = Vrmax * Ec;
  double vrmin_eff = Vrmin * Ec;

  if (zero_TA) {
    double raw = KA * Regulator_blk_in;
    VR = std::min(vrmax_eff, std::max(vrmin_eff, raw));
  } else {
    // Use 8-arg form to pass dynamic limits: (u, dt, xmin, xmax, ymin, ymax, stage, update)
    VR = Regulator_blk.getoutput(Regulator_blk_in, t_inc,
                                  -1000.0, 1000.0,
                                  vrmin_eff, vrmax_eff,
                                  int_flag, true);
  }

  double sat_signal = KE*Efd;
  if (has_Sat) sat_signal += Efd*Sat(Efd);

  double output_blk_in = VR - sat_signal;
  Efd = Output_blk.getoutput(output_blk_in, t_inc, int_flag, true);
}

void gridpack::dynamic_simulation::Ieeet2Model::predictor(
    double t_inc, bool flag)
{
  computeModel(t_inc, PREDICTOR);
}

void gridpack::dynamic_simulation::Ieeet2Model::corrector(
    double t_inc, bool flag)
{
  computeModel(t_inc, CORRECTOR);
}

void gridpack::dynamic_simulation::Ieeet2Model::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

double gridpack::dynamic_simulation::Ieeet2Model::getFieldVoltage()
{
  return Efd;
}

void gridpack::dynamic_simulation::Ieeet2Model::setVterminal(double Vm)
{
  Ec = Vm;
}

void gridpack::dynamic_simulation::Ieeet2Model::setVstab(double Vstab)
{
  Vs = Vstab;
}

void gridpack::dynamic_simulation::Ieeet2Model::setExtBusNum(int ExtBusNum)
{
  p_bus_num = ExtBusNum;
}

void gridpack::dynamic_simulation::Ieeet2Model::setExtGenId(std::string ExtGenId)
{
  p_gen_id = ExtGenId;
}

bool gridpack::dynamic_simulation::Ieeet2Model::setState(
    std::string name, double value)
{
  return false;
}

bool gridpack::dynamic_simulation::Ieeet2Model::getState(
    std::string name, double *value)
{
  return false;
}
