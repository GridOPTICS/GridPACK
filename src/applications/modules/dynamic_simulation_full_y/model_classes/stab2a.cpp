/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   stab2a.cpp
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief  STAB2A — Simple PSS stabilizer (speed-input).
 *
 *  Signal chain:
 *    dw = omega - 1.0
 *    WO = KT*T*s/(1+T*s) * dw        (Cblock: gain KT baked into numerator)
 *    LL1 = (1+T1*s)/(1+T2*s) * WO
 *    LL2 = (1+T3*s)/(1+T4*s) * LL1
 *    Vstab = clamp(LL2, H2, H1)
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cmath>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_pss_model.hpp"
#include "stab2a.hpp"

#define TS_THRESHOLD 4

gridpack::dynamic_simulation::Stab2aModel::Stab2aModel(void)
{
  omega    = 1.0;
  Vstab    = 0.0;
  zero_T2  = false;
  zero_T4  = false;
}

gridpack::dynamic_simulation::Stab2aModel::~Stab2aModel(void)
{
}

void gridpack::dynamic_simulation::Stab2aModel::load(
    boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  if (!data->getValue(STAB2A_KT, &KT, idx)) KT =  1.0;
  if (!data->getValue(STAB2A_T,  &T,  idx)) T  = 10.0;
  if (!data->getValue(STAB2A_T1, &T1, idx)) T1 =  0.0;
  if (!data->getValue(STAB2A_T2, &T2, idx)) T2 =  0.0;
  if (!data->getValue(STAB2A_T3, &T3, idx)) T3 =  0.0;
  if (!data->getValue(STAB2A_T4, &T4, idx)) T4 =  0.0;
  if (!data->getValue(STAB2A_H1, &H1, idx)) H1 =  0.1;
  if (!data->getValue(STAB2A_H2, &H2, idx)) H2 = -0.1;

  if (fabs(T2) < 1e-6) zero_T2 = true;
  if (fabs(T4) < 1e-6) zero_T4 = true;

  // Washout: KT * T*s/(1+T*s)  → Cblock with K=KT*T, T=T
  double wo_a[2], wo_b[2];
  wo_a[0] = T;      wo_a[1] = 1.0;
  wo_b[0] = KT * T; wo_b[1] = 0.0;
  WO_blk.setcoeffs(wo_a, wo_b);

  if (!zero_T2) LL1_blk.setparams(T1, T2);
  if (!zero_T4) LL2_blk.setparams(T3, T4);
}

void gridpack::dynamic_simulation::Stab2aModel::init(
    double mag, double ang, double ts)
{
  double mult_ts = TS_THRESHOLD * ts;

  if (!zero_T2 && T2 > 0.0 && T2 < mult_ts) {
    T2 = mult_ts;
    LL1_blk.setparams(T1, T2);
  }
  if (!zero_T4 && T4 > 0.0 && T4 < mult_ts) {
    T4 = mult_ts;
    LL2_blk.setparams(T3, T4);
  }

  // At steady state: dw=0, all block outputs=0
  if (!zero_T4) LL2_blk.init_given_u(0.0);
  if (!zero_T2) LL1_blk.init_given_u(0.0);
  WO_blk.init_given_u(0.0);

  Vstab = 0.0;
}

void gridpack::dynamic_simulation::Stab2aModel::computeModel(
    double t_inc, IntegrationStage int_flag)
{
  double dw = omega - 1.0;

  // Gain + washout combined: KT * T*s/(1+T*s) * dw
  double wo_out = WO_blk.getoutput(dw, t_inc, int_flag, true);

  // Lead-lag 1: (1+T1*s)/(1+T2*s)
  double ll1_out;
  if (!zero_T2) {
    ll1_out = LL1_blk.getoutput(wo_out, t_inc, int_flag, true);
  } else {
    ll1_out = (fabs(T2) < 1e-6 && fabs(T1) > 1e-6) ? T1 * wo_out : wo_out;
  }

  // Lead-lag 2: (1+T3*s)/(1+T4*s)
  double ll2_out;
  if (!zero_T4) {
    ll2_out = LL2_blk.getoutput(ll1_out, t_inc, int_flag, true);
  } else {
    ll2_out = (fabs(T4) < 1e-6 && fabs(T3) > 1e-6) ? T3 * ll1_out : ll1_out;
  }

  // Output clamp
  if (ll2_out > H1)       Vstab = H1;
  else if (ll2_out < H2)  Vstab = H2;
  else                    Vstab = ll2_out;
}

void gridpack::dynamic_simulation::Stab2aModel::predictor(
    double t_inc, bool flag)
{
  computeModel(t_inc, PREDICTOR);
}

void gridpack::dynamic_simulation::Stab2aModel::corrector(
    double t_inc, bool flag)
{
  computeModel(t_inc, CORRECTOR);
}

double gridpack::dynamic_simulation::Stab2aModel::getVstab()
{
  return Vstab;
}

void gridpack::dynamic_simulation::Stab2aModel::setOmega(double om)
{
  // om is speed deviation (x2w from GENROU); convert to absolute speed
  omega = om + 1.0;
}
