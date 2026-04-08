/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   scrx.cpp
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief  SCRX — Simple Controllable Rectifier Exciter
 *
 *  Same topology as SEXS.  Only difference: when CSWITCH=0 (bus-fed)
 *  the output limits are scaled by the terminal voltage Vt each step.
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_exciter_model.hpp"
#include "scrx.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::ScrxModel::ScrxModel(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::ScrxModel::~ScrxModel(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 */
void gridpack::dynamic_simulation::ScrxModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(EXCITER_TA_OVER_TB, &TA_OVER_TB, idx)) TA_OVER_TB = 0.0;
  if (!data->getValue(EXCITER_TB,         &TB,         idx)) TB         = 0.0;
  if (!data->getValue(EXCITER_K,          &K,          idx)) K          = 0.0;
  if (!data->getValue(EXCITER_TE,         &TE,         idx)) TE         = 0.0;
  if (!data->getValue(EXCITER_EMIN,       &EMIN,       idx)) EMIN       = 0.0;
  if (!data->getValue(EXCITER_EMAX,       &EMAX,       idx)) EMAX       = 999.0;

  // CSWITCH stored under EXCITER_SWITCH: 0=bus-fed (limits×Vt), 1=solid-fed
  double sw_d = 1.0;
  if (!data->getValue(EXCITER_SWITCH, &sw_d, idx)) sw_d = 1.0;
  CSWITCH = (int)sw_d;

  TA = TA_OVER_TB * TB;

  leadlagblock.setparams(TA, TB);

  zero_TE = (fabs(TE) < 1e-6);

  if (!zero_TE) {
    // Limits set here are defaults; they are overridden dynamically when
    // CSWITCH=0, but for CSWITCH=1 these are the effective limits.
    filterblock.setparams(K, TE, EMIN, EMAX, -1000.0, 1000.0);
  } else {
    gainblock.setparams(K, EMIN, EMAX);
  }

  Vs = 0.0;
}

/**
 * Initialize exciter model before calculation
 */
void gridpack::dynamic_simulation::ScrxModel::init(double Vm, double ang, double ts)
{
  double y1;

  // Backward init: start from Efd and work upstream
  // For initialization Ec = Vm; CSWITCH Vt-scaling is Vm ≈ 1.0 pu anyway
  if (!zero_TE) {
    y1 = filterblock.init_given_y(Efd);
  } else {
    y1 = std::min(EMAX, std::max(EMIN, Efd / K));
  }

  double u1 = leadlagblock.init_given_y(y1);

  Vref = Vm + u1 - Vs;
}

/**
 * Predict new state variables for time step
 */
void gridpack::dynamic_simulation::ScrxModel::predictor(double t_inc, bool flag)
{
  double u1 = Vref + Vs - Ec;
  double y1 = leadlagblock.getoutput(u1, t_inc, PREDICTOR, true);

  // Effective output limits: scale by Vt when bus-fed (CSWITCH=0)
  double emin_eff = (CSWITCH == 0) ? EMIN * Ec : EMIN;
  double emax_eff = (CSWITCH == 0) ? EMAX * Ec : EMAX;

  if (!zero_TE) {
    // xmin/xmax wide (don't limit state), ymin/ymax are effective output limits
    Efd = filterblock.getoutput(y1, t_inc, -1000.0, 1000.0, emin_eff, emax_eff, PREDICTOR, true);
  } else {
    Efd = gainblock.getoutput(y1, emin_eff, emax_eff);
  }
}

/**
 * Correct state variables for time step
 */
void gridpack::dynamic_simulation::ScrxModel::corrector(double t_inc, bool flag)
{
  double u1 = Vref + Vs - Ec;
  double y1 = leadlagblock.getoutput(u1, t_inc, CORRECTOR, true);

  double emin_eff = (CSWITCH == 0) ? EMIN * Ec : EMIN;
  double emax_eff = (CSWITCH == 0) ? EMAX * Ec : EMAX;

  if (!zero_TE) {
    Efd = filterblock.getoutput(y1, t_inc, -1000.0, 1000.0, emin_eff, emax_eff, CORRECTOR, true);
  } else {
    Efd = gainblock.getoutput(y1, emin_eff, emax_eff);
  }
}

/**
 * Set the field voltage parameter inside the exciter
 */
void gridpack::dynamic_simulation::ScrxModel::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

/**
 * Get the value of the field voltage parameter
 */
double gridpack::dynamic_simulation::ScrxModel::getFieldVoltage()
{
  return Efd;
}

/**
 * Set the value of terminal voltage
 */
void gridpack::dynamic_simulation::ScrxModel::setVterminal(double Vm)
{
  Ec = Vm;
}

void gridpack::dynamic_simulation::ScrxModel::setVstab(double Vstab)
{
  Vs = Vstab;
}

/**
 * Set the exciter bus number
 */
void gridpack::dynamic_simulation::ScrxModel::setExtBusNum(int ExtBusNum)
{
  p_bus_id = ExtBusNum;
}

/**
 * Set the exciter generator id
 */
void gridpack::dynamic_simulation::ScrxModel::setExtGenId(std::string ExtGenId)
{
  p_ckt = ExtGenId;
}

/**
 * Set internal state parameter in exciter
 */
bool gridpack::dynamic_simulation::ScrxModel::setState(std::string name, double value)
{
  return false;
}

/**
 * Get internal state parameter in exciter
 */
bool gridpack::dynamic_simulation::ScrxModel::getState(std::string name, double *value)
{
  return false;
}
