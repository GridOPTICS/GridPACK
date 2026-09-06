/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   regfma1.hpp
 * @brief  WECC grid-forming converter (PSS/E REGFM_A1) for the EMT module:
 *         droop-controlled voltage source behind Zsource with P/Q limit PIs
 *         and current limiting. Controls are explicit; states are the currents.
 */

#ifndef _regfma1_h_
#define _regfma1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_gen_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

class Regfma1 : public BaseEMTGenModel
{
public:
  Regfma1();
  virtual ~Regfma1();

  void load(boost::shared_ptr<gridpack::component::DataCollection> data, int idx);
  void init(gridpack::RealType *values);
  bool serialWrite(char *string, const int bufsize, const char *signal);
  void write(const char* signal, char* string);
  void getCurrent(double *ia, double *ib, double *ic);
  void getCurrentGlobalLocation(int *i_gloc);
  void getPower(double time, double *Pg, double *Qg);
  void getInitialPower(double *Pg, double *Qg);
  double getFreq();
  double getSpeedDeviation() { return (omega - OMEGA_S)/OMEGA_S; }
  double getSpeedDeviation(int *dw_gloc) { *dw_gloc = -1; return getSpeedDeviation(); }
  double getAngle() { return delta; }
  double getAngle(int *delta_gloc) { *delta_gloc = -1; return delta; }
  void preStep(double time, double timestep);
  void postStep(double time);
  void getnvar(int *nvar);
  int matrixNumValues();
  void matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols);
  void vectorGetValues(gridpack::RealType *values);
  void setValues(gridpack::RealType *values);
  double getInitialFieldVoltage() { return 0.0; }
  double getInitialMechanicalPower() { return 0.0; }

private:
  // Parameters
  double Ra, Xl, L;
  double mq, kpv, kiv, mp, kppmax, kipmax, Pmax, Pmin;
  double Emax, Emin, Tpf, Imax, Qmax, Qmin;
  double kpqmax, kiqmax, Tqf, Tvf;
  int    Vflag, QVflag;
  bool   zero_Tpf, zero_Tqf, zero_Tvf;
  gridpack::ComplexType Zsource;

  // Setpoints from initialization
  double Vset, Pset, Qset;

  // Control blocks and their outputs
  Filter      P_filter_blk, Q_filter_blk, V_filter_blk;
  PIControl   Edroop_PI_blk, Pmax_PI_blk, Pmin_PI_blk, Qmax_PI_blk, Qmin_PI_blk;
  GainLimiter Edroop_limiter_blk;
  Integrator  Delta_blk;
  double Pinv, Qinv, Vmeas, Edroop;
  double Pmax_PI_blk_out, Pmin_PI_blk_out, Qmax_PI_blk_out, Qmin_PI_blk_out;
  double Edroop_max, Edroop_min;
  double omega, delta;

  // States: converter current (machine base) and its system-base copy
  double iabc[3], diabc[3], iout[3];
  double vabc[3], eabc[3], vdq0[3], idq0[3];
};

#endif
