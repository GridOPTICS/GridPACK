/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   exdc2.cpp
 * @author Yousu Chen
 * @Updated:  March 2026
 *
 * @brief  EXDC2 exciter model (IEEE DC2A / PSS/E EXDC2).
 *
 * Identical topology to EXDC1; the voltage amplifier (KA/(1+sTA)) uses
 * anti-windup limits [Vrmin, Vrmax] on the state x, which prevents integrator
 * windup when the output is saturated.  The feedback gain is labeled KF1 in
 * PSS/E EXDC2 (vs KF in EXDC1) but occupies the same DYR parameter slot.
 *
 * PSS/E DYR autocorrections applied (per PowerWorld):
 *   - Swap Vrmax/Vrmin if inverted
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_exciter_model.hpp"
#include "exdc2.hpp"

#define TS_THRESHOLD 4

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Exdc2Model::Exdc2Model(void)
{
  Vs = 0.0;
  sat_A = 0.0;
  sat_B = 0.0;
  has_Sat = true;
  has_leadlag = true;
  zero_TA = false;
  zero_TR = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Exdc2Model::~Exdc2Model(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 */
void gridpack::dynamic_simulation::Exdc2Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(EXCITER_TR,    &TR,    idx)) TR    = 0.0;
  if (!data->getValue(EXCITER_KA,    &KA,    idx)) KA    = 0.0;
  if (!data->getValue(EXCITER_TA,    &TA,    idx)) TA    = 0.0;
  if (!data->getValue(EXCITER_TB,    &TB,    idx)) TB    = 0.0;
  if (!data->getValue(EXCITER_TC,    &TC,    idx)) TC    = 0.0;
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 0.0;
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = 0.0;
  if (!data->getValue(EXCITER_KE,    &KE,    idx)) KE    = 0.0;
  if (!data->getValue(EXCITER_TE,    &TE,    idx)) TE    = 0.0;
  // KF1 is stored under the same key as KF (identical DYR slot)
  if (!data->getValue(EXCITER_KF,    &KF1,   idx)) KF1   = 0.0;
  if (!data->getValue(EXCITER_TF1,   &TF1,   idx)) TF1   = 0.0;
  if (!data->getValue(EXCITER_E1,    &E1,    idx)) E1    = 0.0;
  if (!data->getValue(EXCITER_SE1,   &SE1,   idx)) SE1   = 0.0;
  if (!data->getValue(EXCITER_E2,    &E2,    idx)) E2    = 0.0;
  if (!data->getValue(EXCITER_SE2,   &SE2,   idx)) SE2   = 0.0;

  if (fabs(SE1 * SE2) < 1e-6) has_Sat = false;
  if (fabs(TB)        < 1e-6) has_leadlag = false;

  // Autocorrection: swap Vrmax/Vrmin if inverted
  if (Vrmax < Vrmin) {
    double tmp = Vrmax; Vrmax = Vrmin; Vrmin = tmp;
  }

  if (has_Sat) {
    double r = std::sqrt(SE1 * E1 / (SE2 * E2));
    sat_A = (E1 - r * E2) / (1.0 - r);
    sat_B = SE1 * E1 / ((E1 - sat_A) * (E1 - sat_A));
  }

  // Set up blocks
  if (!zero_TR) {
    Vmeas_blk.setparams(1.0, TR);
  }

  if (has_leadlag) {
    Leadlag_blk.setparams(TC, TB);
  }

  if (!zero_TA) {
    // Anti-windup: state x clamped to [Vrmin, Vrmax]; since y=x for Filter,
    // this prevents integrator windup when output hits its limits.
    Regulator_blk.setparams(KA, TA, Vrmin, Vrmax, -1000.0, 1000.0);
  } else {
    Regulator_gain_blk.setparams(KA, Vrmin, Vrmax);
  }

  double a[2], b[2];
  a[0] = TF1; a[1] = 1.0;
  b[0] = KF1; b[1] = 0.0;
  Feedback_blk.setcoeffs(a, b);

  Output_blk.setparams(TE);
}

/**
 * Saturation function
 * S(x) = B*(x-A)^2/x  if x > A, else 0
 */
double gridpack::dynamic_simulation::Exdc2Model::Sat(double x)
{
  if (x < sat_A) return 0.0;
  return sat_B * (x - sat_A) * (x - sat_A) / x;
}

/**
 * Initialize exciter model before calculation
 */
void gridpack::dynamic_simulation::Exdc2Model::init(double Vm, double Va, double ts)
{
  double mult_ts = TS_THRESHOLD * ts;

  // TR: filter type — zero out if very small
  if (TR > 0.0 && TR < 0.25 * mult_ts) {
    TR = 0.0;
    zero_TR = true;
  } else if (TR > 0.25 * mult_ts && TR < 0.5 * mult_ts) {
    TR = 0.5 * mult_ts;
    Vmeas_blk.setparams(1.0, TR);
  } else if (fabs(TR) < 1e-6) {
    zero_TR = true;
  }

  // TB: lead-lag denominator
  if (TB > 0.0 && TB < 0.5 * mult_ts) {
    TB = 0.0;
    has_leadlag = false;
  } else if (TB > 0.5 * mult_ts && TB < mult_ts) {
    TB = mult_ts;
    if (has_leadlag) Leadlag_blk.setparams(TC, TB);
  }

  // TE
  if (TE > 0.0 && TE < mult_ts) {
    TE = mult_ts;
    Output_blk.setparams(TE);
  }

  // TF1
  if (TF1 > 0.0 && TF1 < mult_ts) {
    TF1 = mult_ts;
    double a[2], b[2];
    a[0] = TF1; a[1] = 1.0;
    b[0] = KF1; b[1] = 0.0;
    Feedback_blk.setcoeffs(a, b);
  }

  // TA
  if (TA > 0.0 && TA < mult_ts) {
    TA = mult_ts;
    zero_TA = false;
    Regulator_blk.setparams(KA, TA, Vrmin, Vrmax, -1000.0, 1000.0);
  } else if (fabs(TA) < 1e-6) {
    zero_TA = true;
  }

  VF = Feedback_blk.init_given_u(Efd);

  double output_blk_in = Output_blk.init_given_y(Efd);

  double sat_signal = KE * Efd;
  if (has_Sat) {
    sat_signal += Efd * Sat(Efd);
  }

  VR = output_blk_in + sat_signal;

  // Widen Vrmax/Vrmin if initial VR is outside bounds
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

/**
 * computeModel - updates all blocks and advances Efd
 */
void gridpack::dynamic_simulation::Exdc2Model::computeModel(
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

  if (zero_TA) {
    VR = Regulator_gain_blk.getoutput(VLL);
  } else {
    VR = Regulator_blk.getoutput(VLL, t_inc, int_flag, true);
  }

  double sat_signal = KE * Efd;
  if (has_Sat) {
    sat_signal += Efd * Sat(Efd);
  }

  double output_blk_in = VR - sat_signal;

  Efd = Output_blk.getoutput(output_blk_in, t_inc, int_flag, true);
}

/**
 * Predict new state variables for time step
 */
void gridpack::dynamic_simulation::Exdc2Model::predictor(double t_inc, bool flag)
{
  computeModel(t_inc, PREDICTOR);
}

/**
 * Correct state variables for time step
 */
void gridpack::dynamic_simulation::Exdc2Model::corrector(double t_inc, bool flag)
{
  computeModel(t_inc, CORRECTOR);
}

/**
 * Set the field voltage parameter inside the exciter
 */
void gridpack::dynamic_simulation::Exdc2Model::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

/**
 * Get the value of the field voltage parameter
 */
double gridpack::dynamic_simulation::Exdc2Model::getFieldVoltage()
{
  return Efd;
}

/**
 * Set the value of terminal voltage
 */
void gridpack::dynamic_simulation::Exdc2Model::setVterminal(double Vm)
{
  Ec = Vm;
}

/**
 * Set stabilizer signal input
 */
void gridpack::dynamic_simulation::Exdc2Model::setVstab(double Vstab)
{
  Vs = Vstab;
}

/**
 * Set the exciter bus number
 */
void gridpack::dynamic_simulation::Exdc2Model::setExtBusNum(int ExtBusNum)
{
  p_bus_num = ExtBusNum;
}

/**
 * Set the exciter generator id
 */
void gridpack::dynamic_simulation::Exdc2Model::setExtGenId(std::string ExtGenId)
{
  p_gen_id = ExtGenId;
}

/**
 * Set internal state parameter in exciter
 */
bool gridpack::dynamic_simulation::Exdc2Model::setState(std::string name,
    double value)
{
  return false;
}

/**
 * Get internal state parameter in exciter
 */
bool gridpack::dynamic_simulation::Exdc2Model::getState(std::string name,
    double *value)
{
  return false;
}
