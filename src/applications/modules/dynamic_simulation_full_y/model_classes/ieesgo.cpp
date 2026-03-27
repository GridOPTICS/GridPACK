/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   ieesgo.cpp
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief  IEESGO steam turbine governor implementation.
 *
 * Model matches ANDES ieesgo.py topology:
 *   F1: Lag K1/(1+T1*s)
 *   F2: LeadLag (1+T2*s)/(1+T3*s)
 *   HL: static clip [PMIN, PMAX]
 *   F3: Lag 1/(1+T4*s)
 *   F4: Lag K2/(1+T5*s)
 *   F5: Lag K3/(1+T6*s)
 *   Pmech = (1-K2)*F3_y + (1-K3)*F4_y + F5_y
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_governor_model.hpp"
#include "ieesgo.hpp"

gridpack::dynamic_simulation::IeesgoModel::IeesgoModel(void)
{
  w     = 0.0;
  Pmech = 1.0;
  Pref0 = 1.0;
}

gridpack::dynamic_simulation::IeesgoModel::~IeesgoModel(void)
{
}

void gridpack::dynamic_simulation::IeesgoModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(GOVERNOR_T1,   &T1,   idx)) T1   = 0.02;
  if (!data->getValue(GOVERNOR_T2,   &T2,   idx)) T2   = 1.0;
  if (!data->getValue(GOVERNOR_T3,   &T3,   idx)) T3   = 1.0;
  if (!data->getValue(GOVERNOR_T4,   &T4,   idx)) T4   = 0.5;
  if (!data->getValue(GOVERNOR_T5,   &T5,   idx)) T5   = 10.0;
  if (!data->getValue(GOVERNOR_T6,   &T6,   idx)) T6   = 0.5;
  if (!data->getValue(GOVERNOR_K1,   &K1,   idx)) K1   = 20.0;
  if (!data->getValue(GOVERNOR_K2,   &K2,   idx)) K2   = 1.0;
  if (!data->getValue(GOVERNOR_K3,   &K3,   idx)) K3   = 1.0;
  if (!data->getValue(GOVERNOR_PMAX, &Pmax, idx)) Pmax = 1.0;
  if (!data->getValue(GOVERNOR_PMIN, &Pmin, idx)) Pmin = 0.0;

  if (Pmax < Pmin) { double tmp = Pmax; Pmax = Pmin; Pmin = tmp; }

  if (T1 != 0) F1_blk.setparams(K1, T1);
  if (T3 != 0) F2_blk.setparams(T2, T3);
  if (T4 != 0) F3_blk.setparams(1.0, T4);
  if (T5 != 0) F4_blk.setparams(K2, T5);
  if (T6 != 0) F5_blk.setparams(K3, T6);
}

void gridpack::dynamic_simulation::IeesgoModel::init(double mag, double ang, double ts)
{
  // Clamp time constants below integration threshold
  double mult_ts = TS_THRESHOLD * ts;
  if (T3 > 0.0 && T3 < 0.25 * mult_ts) {
    T3 = 0.0;
  } else if (T3 > 0.25 * mult_ts && T3 < 0.5 * mult_ts) {
    T3 = 0.5 * mult_ts;
    if (T3 != 0) F2_blk.setparams(T2, T3);
  }

  // At steady state: Pmech = y_HL = Pref0; F1=F2=0
  if (Pmech > Pmax) Pmax = Pmech + 0.1;
  if (Pmech < Pmin) Pmin = Pmech - 0.1;

  Pref0 = Pmech;

  // Backward init of turbine lags
  // F5: y_F5 = K3*K2*Pmech (ss), input = K2*Pmech = y_F4
  double y_F4 = K2 * Pmech;
  double y_F5 = K3 * y_F4;
  if (T6 != 0) F5_blk.init_given_y(y_F5);

  // F4: y_F4 = K2*Pmech (ss), input = Pmech = y_F3
  if (T5 != 0) F4_blk.init_given_y(y_F4);

  // F3: y_F3 = Pmech (ss), input = y_HL = Pmech
  if (T4 != 0) F3_blk.init_given_y(Pmech);

  // F2, F1: at ss speed deviation = 0, so all zero
  if (T3 != 0) F2_blk.init_given_y(0.0);
  if (T1 != 0) F1_blk.init_given_y(0.0);
}

void gridpack::dynamic_simulation::IeesgoModel::computeModel(
    double t_inc, IntegrationStage int_flag)
{
  // F1: lag K1/(1+T1*s)
  double y_F1;
  if (T1 == 0)
    y_F1 = K1 * w;
  else
    y_F1 = F1_blk.getoutput(w, t_inc, int_flag, true);

  // F2: lead-lag (1+T2*s)/(1+T3*s)
  double y_F2;
  if (T3 == 0)
    y_F2 = y_F1;
  else
    y_F2 = F2_blk.getoutput(y_F1, t_inc, int_flag, true);

  // Static limiter: HL = clip(Pref0 - y_F2, Pmin, Pmax)
  double y_HL = Pref0 - y_F2;
  if (y_HL > Pmax) y_HL = Pmax;
  if (y_HL < Pmin) y_HL = Pmin;

  // F3: lag 1/(1+T4*s)
  double y_F3;
  if (T4 == 0)
    y_F3 = y_HL;
  else
    y_F3 = F3_blk.getoutput(y_HL, t_inc, int_flag, true);

  // F4: lag K2/(1+T5*s)
  double y_F4;
  if (T5 == 0)
    y_F4 = K2 * y_F3;
  else
    y_F4 = F4_blk.getoutput(y_F3, t_inc, int_flag, true);

  // F5: lag K3/(1+T6*s)
  double y_F5;
  if (T6 == 0)
    y_F5 = K3 * y_F4;
  else
    y_F5 = F5_blk.getoutput(y_F4, t_inc, int_flag, true);

  // Mechanical power: HP + IP + LP fractions
  Pmech = (1.0 - K2) * y_F3 + (1.0 - K3) * y_F4 + y_F5;
}

void gridpack::dynamic_simulation::IeesgoModel::predictor(double t_inc, bool flag)
{
  computeModel(t_inc, PREDICTOR);
}

void gridpack::dynamic_simulation::IeesgoModel::corrector(double t_inc, bool flag)
{
  computeModel(t_inc, CORRECTOR);
}

void gridpack::dynamic_simulation::IeesgoModel::setMechanicalPower(double pmech)
{
  Pmech = pmech;
}

void gridpack::dynamic_simulation::IeesgoModel::setRotorSpeedDeviation(double dw)
{
  w = dw;
}

double gridpack::dynamic_simulation::IeesgoModel::getMechanicalPower()
{
  return Pmech;
}

bool gridpack::dynamic_simulation::IeesgoModel::setState(std::string name, double value)
{
  return false;
}

bool gridpack::dynamic_simulation::IeesgoModel::getState(std::string name, double *value)
{
  return false;
}
