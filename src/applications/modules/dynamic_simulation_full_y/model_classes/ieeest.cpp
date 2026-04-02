/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   ieeest.cpp
 * @author Yousu Chen
 * @Updated:  March 2026
 *
 * @brief  IEEEST speed-input PSS model.
 *
 * Signal chain:
 *   dw = omega - 1.0  (speed deviation)
 *   F1a = 1/(1+sA1) * dw
 *   F1b = 1/(1+sA2) * F1a
 *   F2a = (1+sA3)/(1+sA4) * F1b
 *   F2b = (1+sA5)/(1+sA6) * F2a
 *   LL1 = (1+sT1)/(1+sT2) * F2b
 *   LL2 = (1+sT3)/(1+sT4) * LL1
 *   WO  = KS*T5*s/(1+sT6) * LL2   [washout, K=KS*T5, T=T6]
 *   Vstab = clamp(WO, LSMIN, LSMAX)
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cmath>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_pss_model.hpp"
#include "ieeest.hpp"

#define TS_THRESHOLD 4

gridpack::dynamic_simulation::IeeeStModel::IeeeStModel(void)
{
  omega = 1.0;
  Vstab = 0.0;
  zero_A1 = false;
  zero_A2 = false;
  zero_A4 = false;
  zero_A6 = false;
  zero_T2 = false;
  zero_T4 = false;
  zero_T6 = false;
}

gridpack::dynamic_simulation::IeeeStModel::~IeeeStModel(void)
{
}

void gridpack::dynamic_simulation::IeeeStModel::load(
    boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  if (!data->getValue(IEEEST_A1,    &A1,    idx)) A1    = 0.0;
  if (!data->getValue(IEEEST_A2,    &A2,    idx)) A2    = 0.0;
  if (!data->getValue(IEEEST_A3,    &A3,    idx)) A3    = 0.0;
  if (!data->getValue(IEEEST_A4,    &A4,    idx)) A4    = 0.0;
  if (!data->getValue(IEEEST_A5,    &A5,    idx)) A5    = 0.0;
  if (!data->getValue(IEEEST_A6,    &A6,    idx)) A6    = 0.0;
  if (!data->getValue(IEEEST_T1,    &T1,    idx)) T1    = 0.0;
  if (!data->getValue(IEEEST_T2,    &T2,    idx)) T2    = 0.0;
  if (!data->getValue(IEEEST_T3,    &T3,    idx)) T3    = 0.0;
  if (!data->getValue(IEEEST_T4,    &T4,    idx)) T4    = 0.0;
  if (!data->getValue(IEEEST_T5,    &T5,    idx)) T5    = 0.0;
  if (!data->getValue(IEEEST_T6,    &T6,    idx)) T6    = 0.0;
  if (!data->getValue(IEEEST_KS,    &KS,    idx)) KS    = 0.0;
  if (!data->getValue(IEEEST_LSMAX, &LSMAX, idx)) LSMAX =  0.1;
  if (!data->getValue(IEEEST_LSMIN, &LSMIN, idx)) LSMIN = -0.1;

  if (fabs(A1) < 1e-6) zero_A1 = true;
  if (fabs(A2) < 1e-6) zero_A2 = true;
  if (fabs(A4) < 1e-6) zero_A4 = true;
  if (fabs(A6) < 1e-6) zero_A6 = true;
  if (fabs(T2) < 1e-6) zero_T2 = true;
  if (fabs(T4) < 1e-6) zero_T4 = true;
  if (fabs(T6) < 1e-6) zero_T6 = true;

  // Set up blocks
  if (!zero_A1) F1a_blk.setparams(1.0, A1);
  if (!zero_A2) F1b_blk.setparams(1.0, A2);
  if (!zero_A4) F2a_blk.setparams(A3, A4);
  if (!zero_A6) F2b_blk.setparams(A5, A6);
  if (!zero_T2) LL1_blk.setparams(T1, T2);
  if (!zero_T4) LL2_blk.setparams(T3, T4);

  // Washout: KS * sT5/(1+sT6) — T5 numerator, T6 denominator (IEEE 421.5 / PW)
  double wo_a[2], wo_b[2];
  wo_a[0] = T6;  wo_a[1] = 1.0;
  wo_b[0] = KS * T5; wo_b[1] = 0.0;
  WO_blk.setcoeffs(wo_a, wo_b);
}

void gridpack::dynamic_simulation::IeeeStModel::init(
    double mag, double ang, double ts)
{
  double mult_ts = TS_THRESHOLD * ts;

  // Time-constant minimum enforcement (same pattern as other models)
  if (!zero_A1 && A1 > 0.0 && A1 < mult_ts) {
    A1 = mult_ts;
    F1a_blk.setparams(1.0, A1);
  }
  if (!zero_A2 && A2 > 0.0 && A2 < mult_ts) {
    A2 = mult_ts;
    F1b_blk.setparams(1.0, A2);
  }
  if (!zero_A4 && A4 > 0.0 && A4 < mult_ts) {
    A4 = mult_ts;
    F2a_blk.setparams(A3, A4);
  }
  if (!zero_A6 && A6 > 0.0 && A6 < mult_ts) {
    A6 = mult_ts;
    F2b_blk.setparams(A5, A6);
  }
  if (!zero_T2 && T2 > 0.0 && T2 < mult_ts) {
    T2 = mult_ts;
    LL1_blk.setparams(T1, T2);
  }
  if (!zero_T4 && T4 > 0.0 && T4 < mult_ts) {
    T4 = mult_ts;
    LL2_blk.setparams(T3, T4);
  }
  if (!zero_T6 && T6 > 0.0 && T6 < mult_ts) {
    T6 = mult_ts;
    // Re-set washout coefficients with corrected T6
    double wo_a2[2], wo_b2[2];
    wo_a2[0] = T6;  wo_a2[1] = 1.0;
    wo_b2[0] = KS * T5; wo_b2[1] = 0.0;
    WO_blk.setcoeffs(wo_a2, wo_b2);
  }

  // At steady state Vstab = 0, so all block outputs are 0.
  // All block states are zero at steady state (speed deviation = 0 at equilibrium).
  WO_blk.init_given_u(0.0);
  if (!zero_T4) LL2_blk.init_given_u(0.0);
  if (!zero_T2) LL1_blk.init_given_u(0.0);
  if (!zero_A6) F2b_blk.init_given_u(0.0);
  if (!zero_A4) F2a_blk.init_given_u(0.0);
  if (!zero_A2) F1b_blk.init_given_u(0.0);
  if (!zero_A1) F1a_blk.init_given_u(0.0);

  Vstab = 0.0;
}

void gridpack::dynamic_simulation::IeeeStModel::computeModel(
    double t_inc, IntegrationStage int_flag)
{
  double dw = omega - 1.0;

  // F1 filter cascade
  if (!zero_A1) {
    F1a = F1a_blk.getoutput(dw,  t_inc, int_flag, true);
  } else {
    F1a = dw;
  }
  if (!zero_A2) {
    F1b = F1b_blk.getoutput(F1a, t_inc, int_flag, true);
  } else {
    F1b = F1a;
  }

  // F2 lead-lag cascade
  if (!zero_A4) {
    F2a = F2a_blk.getoutput(F1b, t_inc, int_flag, true);
  } else {
    F2a = (fabs(A4) < 1e-6 && fabs(A3) > 1e-6) ? A3 * F1b : F1b;
  }
  if (!zero_A6) {
    F2b = F2b_blk.getoutput(F2a, t_inc, int_flag, true);
  } else {
    F2b = (fabs(A6) < 1e-6 && fabs(A5) > 1e-6) ? A5 * F2a : F2a;
  }

  // Lead-lag pair 1
  if (!zero_T2) {
    LL1 = LL1_blk.getoutput(F2b, t_inc, int_flag, true);
  } else {
    LL1 = (fabs(T2) < 1e-6 && fabs(T1) > 1e-6) ? T1 * F2b : F2b;
  }

  // Lead-lag pair 2
  if (!zero_T4) {
    LL2 = LL2_blk.getoutput(LL1, t_inc, int_flag, true);
  } else {
    LL2 = (fabs(T4) < 1e-6 && fabs(T3) > 1e-6) ? T3 * LL1 : LL1;
  }

  // Washout: KS*sT5/(1+sT6) — T5 numerator, T6 denominator
  double WO_out = WO_blk.getoutput(LL2, t_inc, int_flag, true);

  // Clamp to [LSMIN, LSMAX]
  if (WO_out > LSMAX) {
    Vstab = LSMAX;
  } else if (WO_out < LSMIN) {
    Vstab = LSMIN;
  } else {
    Vstab = WO_out;
  }
}

void gridpack::dynamic_simulation::IeeeStModel::predictor(
    double t_inc, bool flag)
{
  computeModel(t_inc, PREDICTOR);
}

void gridpack::dynamic_simulation::IeeeStModel::corrector(
    double t_inc, bool flag)
{
  computeModel(t_inc, CORRECTOR);
}

double gridpack::dynamic_simulation::IeeeStModel::getVstab()
{
  return Vstab;
}

void gridpack::dynamic_simulation::IeeeStModel::setOmega(double om)
{
  // om is speed deviation (x2w from GENROU); convert to absolute speed
  omega = om + 1.0;
}
