/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   gensal.hpp
 * @brief  Salient-pole machine (PSS/E GENSAL) for the EMT module; Genrou
 *         without the q-axis transient state, saturation on E'q.
 */

#ifndef _gensal_model_h_
#define _gensal_model_h_

#include <base_gen_model.hpp>
#include <gridpack/include/gridpack.hpp>

class Gensal: public BaseEMTGenModel
{
public:
  Gensal();
  ~Gensal();

  void load(const boost::shared_ptr<gridpack::component::DataCollection>
            data, int idx);
  bool setJacobian(gridpack::RealType **values);
  void init(gridpack::RealType *values);
  bool serialWrite(char *string, const int bufsize, const char *signal);
  void write(const char* signal, char* string);
  void getCurrent(double *ia, double *ib, double *ic);
  void getCurrentGlobalLocation(int *i_gloc);
  double Sat(double x);
  double dSat(double x);
  double getSpeedDeviation() { return dw; }
  double getSpeedDeviation(int *dw_gloc)
  {
    *dw_gloc = (integrationtype != EXPLICIT) ? p_gloc + 7 : -1;
    return dw;
  }
  void preStep(double time, double timestep);
  void postStep(double time);
  void getnvar(int *nvar);
  int matrixNumValues();
  void matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols);
  void vectorGetValues(gridpack::RealType *values);
  void setValues(gridpack::RealType *values);
  double getInitialFieldVoltage();
  double getFieldCurrent() { return LadIfd; }
  double getInitialMechanicalPower() { return TM; }
  double getAngle() { return delta; }
  double getAngle(int *delta_gloc)
  {
    *delta_gloc = (integrationtype != EXPLICIT) ? p_gloc + 6 : -1;
    return delta;
  }

private:
  double H, D, Ra, Xd, Xq, Xdp, Xdpp, Xl;
  double Tdop, Tdopp, Tqopp, S10, S12;

  double LadIfd;   // Field current
  double TM;       // Mechanical torque
  double VD, VQ;
  double vdq0[3], idq0[3], vabc[3];
  int flux_speed_sensitivity;

  // States: psid psiq psi0 | Eqp psi1d psi2q delta dw | ia ib ic
  double psid, psiq, psi0;
  double Eqp, psi1d, psi2q, delta, dw;
  double iabc[3];
  double dpsid, dpsiq, dpsi0, dEqp, dpsi1d, dpsi2q, ddelta, ddw, diabc[3];

  bool enableSat;
  int bid;
};

#endif
