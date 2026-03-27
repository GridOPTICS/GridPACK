/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   st2cut.cpp
 * @author Yousu Chen
 * @Updated:  March 2026
 *
 * @brief  ST2CUT two-input PSS model.
 */

#include <vector>
#include <iostream>
#include <cmath>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_pss_model.hpp"
#include "st2cut.hpp"

#define TS_THRESHOLD 4

gridpack::dynamic_simulation::St2cutModel::St2cutModel(void)
{
  omega  = 1.0;
  omega2 = 1.0;
  Vstab  = 0.0;
  Vterm  = 1.0;
  VCU    = 999.0;
  VCL    = -999.0;
  zero_T1  = false;
  zero_T2  = false;
  zero_T4  = false;
  zero_T6  = false;
  zero_T8  = false;
  zero_T10 = false;
}

gridpack::dynamic_simulation::St2cutModel::~St2cutModel(void)
{
}

void gridpack::dynamic_simulation::St2cutModel::load(
    boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  if (!data->getValue(ST2CUT_MODE,  &MODE,  idx)) MODE  = 1;
  if (!data->getValue(ST2CUT_MODE2, &MODE2, idx)) MODE2 = 1;
  if (!data->getValue(ST2CUT_K1,    &K1,    idx)) K1    = 0.0;
  if (!data->getValue(ST2CUT_K2,    &K2,    idx)) K2    = 0.0;
  if (!data->getValue(ST2CUT_T1,    &T1,    idx)) T1    = 0.0;
  if (!data->getValue(ST2CUT_T2,    &T2,    idx)) T2    = 0.0;
  if (!data->getValue(ST2CUT_T3,    &T3,    idx)) T3    = 0.0;
  if (!data->getValue(ST2CUT_T4,    &T4,    idx)) T4    = 0.0;
  if (!data->getValue(ST2CUT_T5,    &T5,    idx)) T5    = 0.0;
  if (!data->getValue(ST2CUT_T6,    &T6,    idx)) T6    = 0.0;
  if (!data->getValue(ST2CUT_T7,    &T7,    idx)) T7    = 0.0;
  if (!data->getValue(ST2CUT_T8,    &T8,    idx)) T8    = 0.0;
  if (!data->getValue(ST2CUT_T9,    &T9,    idx)) T9    = 0.0;
  if (!data->getValue(ST2CUT_T10,   &T10,   idx)) T10   = 0.0;
  if (!data->getValue(ST2CUT_LSMAX, &LSMAX, idx)) LSMAX =  0.1;
  if (!data->getValue(ST2CUT_LSMIN, &LSMIN, idx)) LSMIN = -0.1;
  if (!data->getValue(ST2CUT_VCU, &VCU, idx)) VCU = 999.0;
  if (!data->getValue(ST2CUT_VCL, &VCL, idx)) VCL = -999.0;

  if (fabs(T1)  < 1e-6) zero_T1  = true;
  if (fabs(T2)  < 1e-6) zero_T2  = true;
  if (fabs(T4)  < 1e-6) zero_T4  = true;
  if (fabs(T6)  < 1e-6) zero_T6  = true;
  if (fabs(T8)  < 1e-6) zero_T8  = true;
  if (fabs(T10) < 1e-6) zero_T10 = true;

  // Input filter blocks: K/(1+sT)
  if (!zero_T1) In1_blk.setparams(K1, T1);
  if (!zero_T2) In2_blk.setparams(K2, T2);

  // Washout: T3*s/(1+sT4)
  double wo_a[2], wo_b[2];
  wo_a[0] = T4;  wo_a[1] = 1.0;
  wo_b[0] = T3;  wo_b[1] = 0.0;
  WO_blk.setcoeffs(wo_a, wo_b);

  if (!zero_T6)  LL1_blk.setparams(T5, T6);
  if (!zero_T8)  LL2_blk.setparams(T7, T8);
  if (!zero_T10) LL3_blk.setparams(T9, T10);
}

void gridpack::dynamic_simulation::St2cutModel::init(
    double mag, double ang, double ts)
{
  double mult_ts = TS_THRESHOLD * ts;

  if (!zero_T1 && T1 > 0.0 && T1 < mult_ts) {
    T1 = mult_ts;
    In1_blk.setparams(K1, T1);
  }
  if (!zero_T2 && T2 > 0.0 && T2 < mult_ts) {
    T2 = mult_ts;
    In2_blk.setparams(K2, T2);
  }
  if (!zero_T4 && T4 > 0.0 && T4 < mult_ts) {
    T4 = mult_ts;
    double wo_a[2], wo_b[2];
    wo_a[0] = T4;  wo_a[1] = 1.0;
    wo_b[0] = T3;  wo_b[1] = 0.0;
    WO_blk.setcoeffs(wo_a, wo_b);
  }
  if (!zero_T6 && T6 > 0.0 && T6 < mult_ts) {
    T6 = mult_ts;
    LL1_blk.setparams(T5, T6);
  }
  if (!zero_T8 && T8 > 0.0 && T8 < mult_ts) {
    T8 = mult_ts;
    LL2_blk.setparams(T7, T8);
  }
  if (!zero_T10 && T10 > 0.0 && T10 < mult_ts) {
    T10 = mult_ts;
    LL3_blk.setparams(T9, T10);
  }

  // At steady state all inputs are zero (speed deviation = 0)
  WO_blk.init_given_u(0.0);
  if (!zero_T10) LL3_blk.init_given_u(0.0);
  if (!zero_T8)  LL2_blk.init_given_u(0.0);
  if (!zero_T6)  LL1_blk.init_given_u(0.0);
  if (!zero_T2)  In2_blk.init_given_u(0.0);
  if (!zero_T1)  In1_blk.init_given_u(0.0);

  Vstab = 0.0;
}

void gridpack::dynamic_simulation::St2cutModel::computeModel(
    double t_inc, IntegrationStage int_flag)
{
  // Both inputs: speed deviation (pu). MODE selects channel.
  // For now both channels use rotor speed (omega). Only MODE=1 supported.
  double in1 = omega  - 1.0;
  double in2 = omega2 - 1.0;

  // Input filter blocks
  if (!zero_T1) {
    sig1 = In1_blk.getoutput(in1, t_inc, int_flag, true);
  } else {
    sig1 = K1 * in1;
  }
  if (!zero_T2) {
    sig2 = In2_blk.getoutput(in2, t_inc, int_flag, true);
  } else {
    sig2 = K2 * in2;
  }

  double IN = sig1 + sig2;

  // Washout
  double WO_out = WO_blk.getoutput(IN, t_inc, int_flag, true);

  // Lead-lag stages
  if (!zero_T6) {
    LL1 = LL1_blk.getoutput(WO_out, t_inc, int_flag, true);
  } else {
    LL1 = (fabs(T6) < 1e-6 && fabs(T5) > 1e-6) ? T5 * WO_out : WO_out;
  }
  if (!zero_T8) {
    LL2 = LL2_blk.getoutput(LL1, t_inc, int_flag, true);
  } else {
    LL2 = (fabs(T8) < 1e-6 && fabs(T7) > 1e-6) ? T7 * LL1 : LL1;
  }
  if (!zero_T10) {
    LL3 = LL3_blk.getoutput(LL2, t_inc, int_flag, true);
  } else {
    LL3 = (fabs(T10) < 1e-6 && fabs(T9) > 1e-6) ? T9 * LL2 : LL2;
  }

  // Clamp to [LSMIN, LSMAX]
  if (LL3 > LSMAX) {
    Vstab = LSMAX;
  } else if (LL3 < LSMIN) {
    Vstab = LSMIN;
  } else {
    Vstab = LL3;
  }

  // VCU/VCL voltage threshold: disable PSS output when terminal voltage
  // is outside [VCL, VCU] (e.g. during faults with low voltage)
  if (Vterm > VCU || Vterm < VCL) {
    Vstab = 0.0;
  }

}

void gridpack::dynamic_simulation::St2cutModel::predictor(
    double t_inc, bool flag)
{
  computeModel(t_inc, PREDICTOR);
}

void gridpack::dynamic_simulation::St2cutModel::corrector(
    double t_inc, bool flag)
{
  computeModel(t_inc, CORRECTOR);
}

double gridpack::dynamic_simulation::St2cutModel::getVstab()
{
  return Vstab;
}

void gridpack::dynamic_simulation::St2cutModel::setOmega(double om)
{
  // om is speed deviation (x2w from GENROU); convert to absolute speed
  omega  = om + 1.0;
  omega2 = om + 1.0;
}

void gridpack::dynamic_simulation::St2cutModel::setVterminal(double mag)
{
  Vterm = mag;
}
