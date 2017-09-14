/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   ggov1.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 11, 2015
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_governor_model.hpp"
#include "ggov1.hpp"

#define TS_THRESHOLD 1

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Ggov1Model::Ggov1Model(void)
{
  SecondGenExists = false;
  OptionToModifyLimitsForInitialStateLimitViolation = false;
  w = 0.0;
  dx1Pelec = 0;
  dx2GovDer = 0; 
  dx3GovInt = 0;
  dx4Act = 0;
  dx5LL = 0;
  dx6Fload = 0;
  dx7LoadInt = 0;
  dx8LoadCtrl = 0;
  dx9Accel = 0;
  dx10TempLL = 0;
  dx1Pelec_1 = 0;
  dx2GovDer_1 = 0; 
  dx3GovInt_1 = 0;
  dx4Act_1 = 0;
  dx5LL_1 = 0;
  dx6Fload_1 = 0;
  dx7LoadInt_1 = 0;
  dx8LoadCtrl_1 = 0;
  dx9Accel_1 = 0;
  dx10TempLL_1 = 0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Ggov1Model::~Ggov1Model(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to
 * Ggov1Model
 */
void gridpack::dynamic_simulation::Ggov1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  //if (!data->getValue(GOVERNOR_RSELECT, &Rselect, idx)) 
  Rselect = 0.0;
  //if (!data->getValue(GOVERNOR_FLAG, &Flag, idx)) 
  Flag = 0.0;
  //if (!data->getValue(GOVERNOR_R, &R, idx)) 
  R = 0.0; 
  //if (!data->getValue(GOVERNOR_TPELEC, &Tpelec, idx)) 
  Tpelec = 0.0;
  //if (!data->getValue(GOVERNOR_MAXERR, &MaxErr, idx)) 
  MaxErr = 0.0;
  //if (!data->getValue(GOVERNOR_MINERR, &MinErr, idx)) 
  MinErr = 0.0;
  //if (!data->getValue(GOVERNOR_KPGOV, &Kpgov, idx)) 
  Kpgov = 0.0;
  //if (!data->getValue(GOVERNOR_KIGOV, &Kigov, idx)) 
  Kigov = 0.0;
  //if (!data->getValue(GOVERNOR_KDGOV, &Kdgov, idx)) 
  Kdgov = 0.0;
  //if (!data->getValue(GOVERNOR_TDGOV, &Tdgov, idx)) 
  Tdgov = 0.0;
  //if (!data->getValue(GOVERNOR_VMAX, &Vmax, idx)) 
  Vmax = 0.0;
  //if (!data->getValue(GOVERNOR_VMIN, &Vmin, idx)) 
  Vmin = 0.0;
  //if (!data->getValue(GOVERNOR_TACT, &Tact, idx)) 
  Tact = 0.0;
  //if (!data->getValue(GOVERNOR_KTURB, &Kturb, idx)) 
  Kturb = 0.0;
  //if (!data->getValue(GOVERNOR_WFNL, &Wfnl, idx)) 
  Wfnl = 0.0;
  //if (!data->getValue(GOVERNOR_TB, &Tb, idx)) 
  Tb = 0.0; 
  //if (!data->getValue(GOVERNOR_TC, &Tc, idx)) 
  Tc = 0.0; 
  //if (!data->getValue(GOVERNOR_TENG, &Teng, idx)) 
  Teng = 0.0; 
  //if (!data->getValue(GOVERNOR_TFLOAD, &Tfload, idx)) 
  Tfload = 0.0;
  //if (!data->getValue(GOVERNOR_KPLOAD, &Kpload, idx)) 
  Kpload = 0.0;
  //if (!data->getValue(GOVERNOR_KILOAD, &Kiload, idx)) 
  Kiload = 0.0;
  //if (!data->getValue(GOVERNOR_LDREF, &Ldref, idx)) 
  Ldref = 0.0;
  //if (!data->getValue(GOVERNOR_DM, &Dm, idx)) 
  Dm = 0.0;
  //if (!data->getValue(GOVERNOR_ROPEN, &Ropen, idx)) 
  Ropen = 0.0;
  //if (!data->getValue(GOVERNOR_RCLOSE, &Rclose, idx)) 
  Rclose = 0.0;
  //if (!data->getValue(GOVERNOR_KIMW, &Kimw, idx)) 
  Kimw = 0.0;
  //if (!data->getValue(GOVERNOR_ASET, &Aset, idx)) 
  Aset = 0.0;
  //if (!data->getValue(GOVERNOR_KA, &Ka, idx)) 
  Ka = 0.0; 
  //if (!data->getValue(GOVERNOR_TA, &Ta, idx)) 
  Ta = 0.0; 
  //if (!data->getValue(GOVERNOR_TRATE, &Trate, idx)) 
  Trate = 0.0; 
  //if (!data->getValue(GOVERNOR_DB, &Db, idx)) 
  Db = 0.0;
  //if (!data->getValue(GOVERNOR_TSA, &Tsa, idx)) 
  Tsa = 0.0; 
  //if (!data->getValue(GOVERNOR_TSB, &Tsb, idx)) 
  Tsb = 0.0; 
  //if (!data->getValue(GOVERNOR_RUP, &Rup, idx)) 
  Rup = 0.0;
  //if (!data->getValue(GOVERNOR_RDOWN, &Rdown, idx)) 
  Rdown = 0.0;

  if (!data->getValue(GOVERNOR_DB1, &Db1, idx)) Db1 = 0.0; // Db1
  if (!data->getValue(GOVERNOR_ERR, &Err, idx)) Err = 0.0; // Err
  if (!data->getValue(GOVERNOR_DB2, &Db2, idx)) Db2 = 0.0; // Db2

  if (!data->getValue(GENERATOR_MBASE, &GenMVABase, idx)) GenMVABase = 0.0;
}

/**
 * Initialize governor model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Ggov1Model::init(double mag, double ang, double ts)
{
  // Parameter cleanup
  // Just to make the code simpler below we will do the following cleanup 
  // on the input parameters and also store some variables.
  if (Tfload < TS_THRESHOLD * ts) Tfload = TS_THRESHOLD * ts; // force non-zero
  if (Tdgov == 0) Kdgov = 0; // disable derivative
  // Following is done because of deadband. It is less confusing to implement a deadband
  // that starts with a zero input. But the integrator is so slow it never moves.
  if (Kigov == 0) Kigov = 0.00001;
  if (Kpgov != 0) KigovKpgov = Kigov / Kpgov;
  else KigovKpgov = 0;
  if (Tb == 0) Tc = 0; // Disable Lead Lag
  if (Kiload != 0) KiLoadKpLoad = Kiload / Kpload;
  else KiLoadKpLoad = 0;
  if (Kturb != 0) LdRefslashKturb = Ldref / Kturb; // TBD: fKturb; // Ldref/Kturb ?
  else LdRefslashKturb = 0;
  if (Ta <= 0) Ka = 0; // disable rate limiter
  if (Trate == 0) Trate = GenMVABase;
  if (Db < 0) Db = 0; // disable deadband
  if (Teng < 0) Teng = 0; // Actually for now we ignore this here anyway

  printf("ggov1: Pmech = %f\n", Pmech);
  // State 1
  x1Pelec = GenPelec * GenMVABase /Trate;
  Pmwset = x1Pelec;
  // State 5
  double PmechBase = Pmech * GenMVABase / Trate;
  if (Tb < TS_THRESHOLD * ts) x5LL = 0;
  else x5LL = PmechBase *(1 - Tc / Tb);
  // Note we are ignoring the engine delay Teng
  // State 4
  x4Act = PmechBase / Kturb + Wfnl;
  // Artifically increase Kturb if x4Act violated Vmax
  if (x4Act > Vmax && x4Act > Wfnl) {
    double Kturbnew = PmechBase / (Vmax - Wfnl);
    LdRefslashKturb = LdRefslashKturb * Kturb / Kturbnew;
    Kturb = Kturbnew;
    x4Act = Vmax;
  }
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (x4Act > Vmax) Vmax = x4Act; // really handled by change to Kturb
    if (x4Act < Vmin) Vmin = x4Act; 
  }
  // State 10
  if (Tsb < TS_THRESHOLD * ts) x10TempLL = 0;
  else x10TempLL = x4Act *(1 - Tsa / Tsb);
  // State 6
  x6Fload = x4Act;
  // State 7
  x7LoadInt = x4Act;
  // State 9
  x9Accel = 0;
  // Pref
  // Note: Because we force Kigov > 0 we assume input to PID is zero
  if (Rselect == 1) Pref = R * x1Pelec;
  else if (Rselect < 0) Pref = R * x4Act;
  else Pref = 0;
  // State 2
  x2GovDer = 0;
  // State 3
  x3GovInt = 0;
  // State 8
  x8LoadCtrl = 0;
  LastLowValueSelect = x4Act;
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Ggov1Model::predictor(double t_inc, bool flag)
{
  if (!flag) {
    x1Pelec = x1Pelec_1;
    x2GovDer = x2GovDer_1; 
    x3GovInt = x3GovInt_1;
    x4Act = x4Act_1;
    x5LL = x5LL_1;
    x6Fload = x6Fload_1;
    x7LoadInt = x7LoadInt_1;
    x8LoadCtrl = x8LoadCtrl_1;
    x9Accel = x9Accel_1;
    x10TempLL = x10TempLL_1;
  }
  // State 1
  if (Tpelec < TS_THRESHOLD * t_inc) {
    x1Pelec = GenPelec * GenMVABase / Trate;
    dx1Pelec = 0;
  } else 
    dx1Pelec = (GenPelec * GenMVABase / Trate - x1Pelec) / Tpelec;
  // State 8
  if (x8LoadCtrl > +1.1 * R) x8LoadCtrl = +1.1 * R;
  else if (x8LoadCtrl < -1.1 * R) x8LoadCtrl = -1.1 * R;
  dx8LoadCtrl = Kimw * (Pmwset - x1Pelec); // TBD: Kimw(Pmwset - x1Pelec)?
  if (dx8LoadCtrl > 0 && x8LoadCtrl >= 1.1 * R) dx8LoadCtrl = 0;
  if (dx8LoadCtrl < 0 && x8LoadCtrl <= -1.1 * R) dx8LoadCtrl = 0;
  // State 9
  if (Ka > 0) dx9Accel = w * Ta + x9Accel;
  else dx9Accel = 0;
  // State 10 and State 6
  // Note: if Tact = 0, then we are really using the previous timestep
  // value of x4Act at this point
  double TempIn = x4Act;
  if (Flag == 1) TempIn = TempIn * (1 + w);
  if (Dm < 0) TempIn = TempIn * pow(1+w, Dm); 
  if (Tsb < TS_THRESHOLD * t_inc) {
    dx10TempLL = 0;
    dx6Fload = (TempIn - x6Fload) / Tfload; // We require Tfload > 0
  } else {
    dx10TempLL = (TempIn * (1 - Tsa / Tsb) - x10TempLL) / Tsb;
    dx6Fload = (x10TempLL - x6Fload) / Tfload;
  }
  // State 7
  // Note: If Kpload = 0, then integrator [Kiload/s] is fed by
  // input to Kpload and output toward fsrn is used
  if (Kpload != 0)
    dx7LoadInt = KiLoadKpLoad * (LastLowValueSelect - x7LoadInt);
  else
    dx7LoadInt = Kiload * (LdRefslashKturb + Wfnl - x6Fload);
  // Find the input to PID controller
  // Note: feedback of x4Act and LastLowValueSelect is really
  // the previous time step here
  double PIDIn = Pref + x8LoadCtrl - w - PIDIn;
  if (Rselect == +1) PIDIn = PIDIn - R * x1Pelec;
  else if (Rselect == -1) PIDIn = PIDIn - R * x4Act;
  else if (Rselect == -2) PIDIn = PIDIn - R * LastLowValueSelect;
  if (PIDIn > +Db) PIDIn = PIDIn - Db;
  else if (PIDIn < -Db) PIDIn = PIDIn + Db;
  else PIDIn = 0;
  if (PIDIn > MaxErr) PIDIn = MaxErr;
  else if (PIDIn < MinErr) PIDIn = MinErr;
  // State 3
  // Note: If Kpgov = 0, then the integrator block moves back
  // in parallel with the derivative and proportional control 
  if (Kpgov != 0)
    dx3GovInt = KigovKpgov * (LastLowValueSelect - x3GovInt);
  else
    dx3GovInt = Kigov * PIDIn;
  // State 2
  double DerivOut;
  if (Tdgov > 0) {
    dx2GovDer = (PIDIn * (-Kdgov / Tdgov) - x2GovDer) / Tdgov;
    DerivOut = x2GovDer + PIDIn * Kdgov / Tdgov;
  } else {
    dx2GovDer = 0;
    DerivOut = 0;
  }
  // Get New Low Value Select
  double Fsrt = x7LoadInt + Kpload * (LdRefslashKturb + Wfnl - x6Fload); 
  if (Fsrt > 1.0) Fsrt = 1.0;
  double Fsra;
  if (Ka > 0)
    Fsra = (Aset - x9Accel) * Ka * t_inc + LastLowValueSelect; // TBD: t_inc is TimeStep?
  else
    Fsra = 999999;
  // Finally, now update the LastLowValueSelect.
  // Note everything up until this point was using last time step value
  LastLowValueSelect = DerivOut + x3GovInt + PIDIn * Kpgov;
  if (LastLowValueSelect > Fsra) LastLowValueSelect = Fsra;
  if (LastLowValueSelect > Fsrt) LastLowValueSelect = Fsrt;
  if (LastLowValueSelect > Vmax) LastLowValueSelect = Vmax;
  if (LastLowValueSelect < Vmin) LastLowValueSelect = Vmin;
  // State 4
  if (Tact < TS_THRESHOLD * t_inc) {
    dx4Act = 0;
    x4Act = LastLowValueSelect;
  } else {
    dx4Act = (LastLowValueSelect - x4Act) / Tact;
    if (dx4Act > Ropen) dx4Act = Ropen;
    if (dx4Act < Rclose) dx4Act = Rclose;
  }
  // State 5
  TempIn = x4Act;
  if (Flag == 1) TempIn = TempIn * (1 + w);
  TempIn = (TempIn - Wfnl) * Kturb;
  // Note: We are ignoring the Engine Delay here
  if (Tb > TS_THRESHOLD) {
    dx5LL = 0;
    LeadLagOut = TempIn;
  } else {
    dx5LL = (TempIn * (1 - Ta / Tb) - x5LL) / Tb;
    LeadLagOut = TempIn * Ta / Tb + x5LL;
  }

  x1Pelec_1 = x1Pelec + dx1Pelec * t_inc;
  x2GovDer_1 = x2GovDer + dx2GovDer * t_inc;
  x3GovInt_1 = x3GovInt + dx3GovInt * t_inc;
  x4Act_1 = x4Act + dx4Act * t_inc;
  x5LL_1 = x5LL + dx5LL * t_inc;
  x6Fload_1 = x6Fload + dx6Fload * t_inc;
  x7LoadInt_1 = x7LoadInt + dx7LoadInt * t_inc;
  x8LoadCtrl_1 = x8LoadCtrl + dx8LoadCtrl * t_inc;
  x9Accel_1 = x9Accel + dx9Accel * t_inc;
  x10TempLL_1 = x10TempLL + dx10TempLL * t_inc;

  printf("ggov1 dx: %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", dx1Pelec, dx2GovDer, dx3GovInt, dx4Act, dx5LL, dx6Fload, dx7LoadInt, dx8LoadCtrl, dx9Accel, dx10TempLL);
  printf("ggov1 x: %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1Pelec_1, x2GovDer_1, x3GovInt_1, x4Act_1, x5LL_1, x6Fload_1, x7LoadInt_1, x8LoadCtrl_1, x9Accel_1, x10TempLL_1);

  if (Dm > 0) {
    Pmech = (LeadLagOut - Dm * w) * Trate / GenMVABase;
  } else {
    Pmech = LeadLagOut * Trate / GenMVABase;
  } 
  
  printf("ggov1 Pmech = %f\n", Pmech);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Ggov1Model::corrector(double t_inc, bool flag)
{
  // State 1
  if (Tpelec < TS_THRESHOLD * t_inc) {
    x1Pelec_1 = GenPelec * GenMVABase / Trate;
    dx1Pelec_1 = 0;
  } else 
    dx1Pelec_1 = (GenPelec * GenMVABase / Trate - x1Pelec_1) / Tpelec;
  // State 8
  if (x8LoadCtrl_1 > +1.1 * R) x8LoadCtrl_1 = +1.1 * R;
  else if (x8LoadCtrl_1 < -1.1 * R) x8LoadCtrl_1 = -1.1 * R;
  dx8LoadCtrl_1 = Kimw * (Pmwset - x1Pelec_1); // TBD: Kimw(Pmwset - x1Pelec)?
  if (dx8LoadCtrl_1 > 0 && x8LoadCtrl_1 >= 1.1 * R) dx8LoadCtrl_1 = 0;
  if (dx8LoadCtrl_1 < 0 && x8LoadCtrl_1 <= -1.1 * R) dx8LoadCtrl_1 = 0;
  // State 9
  if (Ka > 0) dx9Accel_1 = w * Ta + x9Accel_1;
  else dx9Accel_1 = 0;
  // State 10 and State 6
  // Note: if Tact = 0, then we are really using the previous timestep
  // value of x4Act at this point
  double TempIn = x4Act_1;
  if (Flag == 1) TempIn = TempIn * (1 + w);
  if (Dm < 0) TempIn = TempIn * pow(1+w, Dm); 
  if (Tsb < TS_THRESHOLD * t_inc) {
    dx10TempLL_1 = 0;
    dx6Fload_1 = (TempIn - x6Fload_1) / Tfload; // We require Tfload > 0
  } else {
    dx10TempLL_1 = (TempIn * (1 - Tsa / Tsb) - x10TempLL_1) / Tsb;
    dx6Fload_1 = (x10TempLL_1 - x6Fload_1) / Tfload;
  }
  // State 7
  // Note: If Kpload = 0, then integrator [Kiload/s] is fed by
  // input to Kpload and output toward fsrn is used
  if (Kpload != 0)
    dx7LoadInt_1 = KiLoadKpLoad * (LastLowValueSelect - x7LoadInt_1);
  else
    dx7LoadInt_1 = Kiload * (LdRefslashKturb + Wfnl - x6Fload_1);
  // Find the input to PID controller
  // Note: feedback of x4Act and LastLowValueSelect is really
  // the previous time step here
  double PIDIn = Pref + x8LoadCtrl_1 - w - PIDIn;
  if (Rselect == +1) PIDIn = PIDIn - R * x1Pelec_1;
  else if (Rselect == -1) PIDIn = PIDIn - R * x4Act_1;
  else if (Rselect == -2) PIDIn = PIDIn - R * LastLowValueSelect;
  if (PIDIn > +Db) PIDIn = PIDIn - Db;
  else if (PIDIn < -Db) PIDIn = PIDIn + Db;
  else PIDIn = 0;
  if (PIDIn > MaxErr) PIDIn = MaxErr;
  else if (PIDIn < MinErr) PIDIn = MinErr;
  // State 3
  // Note: If Kpgov = 0, then the integrator block moves back
  // in parallel with the derivative and proportional control 
  if (Kpgov != 0)
    dx3GovInt_1 = KigovKpgov * (LastLowValueSelect - x3GovInt_1);
  else
    dx3GovInt_1 = Kigov * PIDIn;
  // State 2
  double DerivOut;
  if (Tdgov > 0) {
    dx2GovDer_1 = (PIDIn * (-Kdgov / Tdgov) - x2GovDer_1) / Tdgov;
    DerivOut = x2GovDer_1 + PIDIn * Kdgov / Tdgov;
  } else {
    dx2GovDer_1 = 0;
    DerivOut = 0;
  }
  // Get New Low Value Select
  double Fsrt = x7LoadInt_1 + Kpload * (LdRefslashKturb + Wfnl - x6Fload_1); 
  if (Fsrt > 1.0) Fsrt = 1.0;
  double Fsra;
  if (Ka > 0)
    Fsra = (Aset - x9Accel_1) * Ka * t_inc + LastLowValueSelect; // TBD: t_inc is TimeStep?
  else
    Fsra = 999999;
  // Finally, now update the LastLowValueSelect.
  // Note everything up until this point was using last time step value
  LastLowValueSelect = DerivOut + x3GovInt_1 + PIDIn * Kpgov;
  if (LastLowValueSelect > Fsra) LastLowValueSelect = Fsra;
  if (LastLowValueSelect > Fsrt) LastLowValueSelect = Fsrt;
  if (LastLowValueSelect > Vmax) LastLowValueSelect = Vmax;
  if (LastLowValueSelect < Vmin) LastLowValueSelect = Vmin;
  // State 4
  if (Tact < TS_THRESHOLD * t_inc) {
    dx4Act_1 = 0;
    x4Act_1 = LastLowValueSelect;
  } else {
    dx4Act_1 = (LastLowValueSelect - x4Act_1) / Tact;
    if (dx4Act_1 > Ropen) dx4Act_1 = Ropen;
    if (dx4Act_1 < Rclose) dx4Act = Rclose;
  }
  // State 5
  TempIn = x4Act_1;
  if (Flag == 1) TempIn = TempIn * (1 + w);
  TempIn = (TempIn - Wfnl) * Kturb;
  // Note: We are ignoring the Engine Delay here
  if (Tb > TS_THRESHOLD) {
    dx5LL_1 = 0;
    LeadLagOut = TempIn;
  } else {
    dx5LL_1 = (TempIn * (1 - Ta / Tb) - x5LL_1) / Tb;
    LeadLagOut = TempIn * Ta / Tb + x5LL_1;
  }

  x1Pelec_1 = x1Pelec + (dx1Pelec + dx1Pelec_1) / 2.0 * t_inc;
  x2GovDer_1 = x2GovDer + (dx2GovDer + dx2GovDer_1) / 2.0 * t_inc;
  x3GovInt_1 = x3GovInt + (dx3GovInt + dx3GovInt_1) / 2.0 * t_inc;
  x4Act_1 = x4Act + (dx4Act + dx4Act_1) / 2.0 * t_inc;
  x5LL_1 = x5LL + (dx5LL + dx5LL_1) / 2.0 * t_inc;
  x6Fload_1 = x6Fload + (dx6Fload + dx6Fload_1) / 2.0 * t_inc;
  x7LoadInt_1 = x7LoadInt + (dx7LoadInt + dx7LoadInt_1) / 2.0 * t_inc;
  x8LoadCtrl_1 = x8LoadCtrl + (dx8LoadCtrl + dx8LoadCtrl_1) / 2.0 * t_inc;
  x9Accel_1 = x9Accel + (dx9Accel + dx9Accel_1) / 2.0 * t_inc;
  x10TempLL_1 = x10TempLL + (dx10TempLL + dx10TempLL_1) / 2.0 * t_inc;

  printf("ggov1 dx: %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", dx1Pelec, dx2GovDer, dx3GovInt, dx4Act, dx5LL, dx6Fload, dx7LoadInt, dx8LoadCtrl, dx9Accel, dx10TempLL);
  printf("ggov1 x: %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1Pelec_1, x2GovDer_1, x3GovInt_1, x4Act_1, x5LL_1, x6Fload_1, x7LoadInt_1, x8LoadCtrl_1, x9Accel_1, x10TempLL_1);

  if (Dm > 0) {
    Pmech = (LeadLagOut - Dm * w) * Trate / GenMVABase;
  } else {
    Pmech = LeadLagOut * Trate / GenMVABase;
  } 
  
  printf("ggov1 Pmech = %f\n", Pmech);
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power
 */
void gridpack::dynamic_simulation::Ggov1Model::setMechanicalPower(double pmech)
{
  Pmech = pmech; 
}

/**
 * Set the rotor speed deviation inside the governor
 * @param delta_o value of the rotor speed deviation
 */
void gridpack::dynamic_simulation::Ggov1Model::setRotorSpeedDeviation(double delta_o)
{
  w = delta_o;
}

/** 
 * Get the value of the mechanical power
 * @return value of mechanical power
 */
double gridpack::dynamic_simulation::Ggov1Model::getMechanicalPower()
{
  return Pmech; 
}

/** 
 * Get the value of the rotor speed deviation
 * 
 */
/*double gridpack::dynamic_simulation::Ggov1Model::getRotorSpeedDeviation()
{
  return w;
}*/
