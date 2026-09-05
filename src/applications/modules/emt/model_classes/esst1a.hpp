/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   esst1a.hpp
 * @brief  IEEE ST1A static exciter (PSS/E ESST1A) for the EMT module.
 *         States: Vmeas xLL1 xLL2 Va xf. UEL/OEL inputs are not modeled.
 */

#ifndef _esst1a_model_h_
#define _esst1a_model_h_

#include <base_exc_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "cblock.hpp"
#include "dblock.hpp"

class Esst1a: public BaseEMTExcModel
{
public:
  Esst1a();
  ~Esst1a();

  void load(const boost::shared_ptr<gridpack::component::DataCollection>
            data, int idx);
  void init(gridpack::RealType *values);
  bool serialWrite(char *string, const int bufsize, const char *signal);
  void write(const char* signal, char* string);
  void getnvar(int *nvar);
  void preStep(double time, double timestep);
  void postStep(double time);
  int matrixNumValues();
  void matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols);
  void vectorGetValues(gridpack::RealType *values);
  void setValues(gridpack::RealType *values);
  double getFieldVoltage();
  double getFieldVoltage(int *Efd_gloc);
  bool getFieldVoltagePartialDerivatives(int *xexc_loc, double *dEfd_dxexc, double *dEfd_dxgen);
  void setEvent(gridpack::math::RealDAESolver::EventManagerPtr);
  void resetEventFlags();
  void eventFunction(const double&t, gridpack::RealType *state, std::vector<gridpack::RealType >& evalues);
  void eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state);

private:
  // Field voltage from regulator output, terminal voltage and field current.
  double computeEfd(double Va_in, double Vt, double Ifd, bool *clamped);
  double terminalVoltage();

  int    UEL, VOS;
  double TR, VImax, VImin, TC, TB, TC1, TB1;
  double KA, TA, VAmax, VAmin, VRmax, VRmin;
  double KC, KF, TF, KLR, ILR;

  double Vmeas, xLL1, xLL2, Va, xf;
  double dVmeas, dxLL1, dxLL2, dVa, dxf;
  double Efd, Ec, Vref, Vs, VF, Ifd;

  bool zero_TR, zero_TB, zero_TB1, zero_TA, zero_TF;
  bool Vi_at_min, Vi_at_max, Va_at_min, Va_at_max;

  Filter      Vmeas_blk;
  LeadLag     LL1_blk, LL2_blk;
  Filter      Regulator_blk;
  GainLimiter Regulator_gain_blk;
  Cblock      Feedback_blk;
};

class Esst1aEvent
  : public gridpack::math::RealDAESolver::Event
{
public:
  Esst1aEvent(Esst1a *exc): gridpack::math::RealDAESolver::Event(4), p_exc(exc)
  {
    std::fill(p_term.begin(), p_term.end(), false);
    std::fill(p_dir.begin(), p_dir.end(), gridpack::math::CrossZeroNegative);
  }
  ~Esst1aEvent(void) {}
protected:
  Esst1a *p_exc;
  void p_update(const double& t, gridpack::RealType *state);
  void p_handle(const bool *triggered, const double& t, gridpack::RealType *state);
};

#endif
