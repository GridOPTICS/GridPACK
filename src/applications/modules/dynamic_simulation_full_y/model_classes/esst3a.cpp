/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   esst3a.cpp
 *
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief  ESST3A static exciter model (IEEE Std 421.5-2016 Type ST3A).
 *
 * Signal flow:
 *   Vmeas  = 1/(1+sTR) * Vt
 *   Verr   = Vref + Vs - Vmeas
 *   VI     = clamp(Verr, VIMIN, VIMAX)
 *   VLL    = (1+sTC)/(1+sTB) * VI            [lead-lag; bypassed if TB=0]
 *   VR     = KA/(1+sTA) * VLL                [anti-windup VRMIN/VRMAX]
 *   VG     = clamp(KG * Efd, 0, VGMAX)
 *   vrs    = VR - VG
 *   VM     = KM/(1+sTM) * vrs                [anti-windup VMMIN/VMMAX]
 *   VE     = |KPC*(Vd+jVq) + j*(KI+KPC_re*XL)*(Id+jIq)|
 *            where KPC = KP * exp(j*THETAP_rad)
 *   IN     = KC*LadIfd/VE
 *   FEX    = piecewise commutation factor
 *   VB     = clamp(VE*FEX, 0, VBMAX)
 *   Efd    = VM * VB
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <cmath>
#include <algorithm>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_exciter_model.hpp"
#include "esst3a.hpp"

#define TS_THRESHOLD 4

gridpack::dynamic_simulation::Esst3aModel::Esst3aModel(void)
{
  Vs = 0.0;
  LadIfd = 0.0;
  Vd = 0.0; Vq = 1.0;
  Id = 0.0; Iq = 0.0;
  zero_TR = false;
  zero_TB = false;
  zero_TA = false;
  zero_TM = false;
}

gridpack::dynamic_simulation::Esst3aModel::~Esst3aModel(void)
{
}

void gridpack::dynamic_simulation::Esst3aModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(EXCITER_TR,    &TR,    idx)) TR    = 0.0;
  if (!data->getValue(EXCITER_VIMAX, &VIMAX, idx)) VIMAX = 1.0;
  if (!data->getValue(EXCITER_VIMIN, &VIMIN, idx)) VIMIN = -1.0;
  if (!data->getValue(EXCITER_KM,    &KM,    idx)) KM    = 1.0;
  if (!data->getValue(EXCITER_TC,    &TC,    idx)) TC    = 0.0;
  if (!data->getValue(EXCITER_TB,    &TB,    idx)) TB    = 0.0;
  if (!data->getValue(EXCITER_KA,    &KA,    idx)) KA    = 1.0;
  if (!data->getValue(EXCITER_TA,    &TA,    idx)) TA    = 0.0;
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 99.0;
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = -99.0;
  if (!data->getValue(EXCITER_KG,    &KG,    idx)) KG    = 1.0;
  if (!data->getValue(EXCITER_KP,    &KP,    idx)) KP    = 1.0;
  if (!data->getValue(EXCITER_KI,    &KI,    idx)) KI    = 0.0;
  if (!data->getValue(EXCITER_VBMAX, &VBMAX, idx)) VBMAX = 10.0;
  if (!data->getValue(EXCITER_KC,    &KC,    idx)) KC    = 0.0;
  if (!data->getValue(EXCITER_XL,    &XL,    idx)) XL    = 0.0;
  if (!data->getValue(EXCITER_VGMAX, &VGMAX, idx)) VGMAX = 10.0;
  if (!data->getValue(EXCITER_THETAP,&THETAP,idx)) THETAP = 0.0;
  if (!data->getValue(EXCITER_TM,    &TM,    idx)) TM    = 0.0;
  if (!data->getValue(EXCITER_VMMAX, &Vmmax, idx)) Vmmax = 99.0;
  if (!data->getValue(EXCITER_VMMIN, &Vmmin, idx)) Vmmin = 0.0;

  // Autocorrection: swap limits if inverted
  if (Vrmax < Vrmin) { double tmp = Vrmax; Vrmax = Vrmin; Vrmin = tmp; }
  if (VIMAX < VIMIN) { double tmp = VIMAX; VIMAX = VIMIN; VIMIN = tmp; }
  if (Vmmax < Vmmin) { double tmp = Vmmax; Vmmax = Vmmin; Vmmin = tmp; }

  if (fabs(TR) < 1e-6) zero_TR = true;
  if (fabs(TB) < 1e-6) zero_TB = true;
  if (fabs(TA) < 1e-6) zero_TA = true;
  if (fabs(TM) < 1e-6) zero_TM = true;

  // KPC = KP * exp(j * THETAP_rad)
  THETAP_rad = THETAP * M_PI / 180.0;
  KPC_re = KP * cos(THETAP_rad);
  KPC_im = KP * sin(THETAP_rad);

  // Set up blocks
  if (!zero_TR) {
    Vmeas_blk.setparams(1.0, TR);
  }
  if (!zero_TB) {
    Leadlag_blk.setparams(TC, TB);
  }
  if (!zero_TA) {
    Regulator_blk.setparams(KA, TA, Vrmin, Vrmax, -1000.0, 1000.0);
  } else {
    Regulator_gain_blk.setparams(KA, Vrmin, Vrmax);
  }
  if (!zero_TM) {
    InnerReg_blk.setparams(KM, TM, Vmmin, Vmmax, -1000.0, 1000.0);
  } else {
    InnerReg_gain_blk.setparams(KM, Vmmin, Vmmax);
  }
}

// Piecewise commutation factor FEX (per IEEE 421.5)
double gridpack::dynamic_simulation::Esst3aModel::computeFex(double IN) const
{
  if (IN <= 0.0) {
    return 1.0;
  } else if (IN < 0.433) {
    return 1.0 - 0.577 * IN;
  } else if (IN < 0.75) {
    double val = 0.75 - IN * IN;
    return (val > 0.0) ? sqrt(val) : 0.0;
  } else if (IN < 1.0) {
    return 1.732 * (1.0 - IN);
  } else {
    return 0.0;
  }
}

void gridpack::dynamic_simulation::Esst3aModel::init(double Vm, double Va, double ts)
{
  double mult_ts = TS_THRESHOLD * ts;

  // TR
  if (TR > 0.0 && TR < 0.25 * mult_ts) {
    TR = 0.0; zero_TR = true;
  } else if (TR > 0.25 * mult_ts && TR < 0.5 * mult_ts) {
    TR = 0.5 * mult_ts;
    Vmeas_blk.setparams(1.0, TR);
  } else if (fabs(TR) < 1e-6) {
    zero_TR = true;
  }

  // TB
  if (TB > 0.0 && TB < 0.5 * mult_ts) {
    TB = 0.0; zero_TB = true;
  } else if (TB > 0.5 * mult_ts && TB < mult_ts) {
    TB = mult_ts;
    if (!zero_TB) Leadlag_blk.setparams(TC, TB);
  }

  // TA
  if (TA > 0.0 && TA < mult_ts) {
    TA = mult_ts; zero_TA = false;
    Regulator_blk.setparams(KA, TA, Vrmin, Vrmax, -1000.0, 1000.0);
  } else if (fabs(TA) < 1e-6) {
    zero_TA = true;
  }

  // TM
  if (TM > 0.0 && TM < mult_ts) {
    TM = mult_ts; zero_TM = false;
    InnerReg_blk.setparams(KM, TM, Vmmin, Vmmax, -1000.0, 1000.0);
  } else if (fabs(TM) < 1e-6) {
    zero_TM = true;
  }

  // Compute VE from d-q quantities
  // KPC*(Vd+jVq): re = KPC_re*Vd - KPC_im*Vq, im = KPC_im*Vd + KPC_re*Vq
  double ve_re = KPC_re * Vd - KPC_im * Vq + (KI + KPC_im * XL) * (-Iq) - KPC_re * XL * Id;
  // Wait, let me re-derive:
  // VE = KPC*(Vd+jVq) + j*(KI+KPC_re*XL)*(Id+jIq)
  //    = (KPC_re*Vd - KPC_im*Vq) + j*(KPC_im*Vd + KPC_re*Vq)
  //      + j*(KI + KPC_re*XL)*Id + j^2*(KI + KPC_re*XL)*Iq
  // But note: KPC*XL includes KPC_im*XL too in the imaginary coefficient
  // Full formula from IEEE 421.5: VE = KPC*(Vd+jVq) + j*(KI + KPC*XL)*(Id+jIq)
  // where KPC*XL = (KPC_re + j*KPC_im)*XL
  // j*(KI + (KPC_re+j*KPC_im)*XL) = j*(KI + KPC_re*XL) + j^2*KPC_im*XL
  //                                 = j*(KI + KPC_re*XL) - KPC_im*XL
  // Full VE (real, imag):
  // re = KPC_re*Vd - KPC_im*Vq   + (-(KI+KPC_re*XL)*Iq) + (-KPC_im*XL)*Id
  //    = KPC_re*Vd - KPC_im*Vq - (KI+KPC_re*XL)*Iq - KPC_im*XL*Id
  // im = KPC_im*Vd + KPC_re*Vq   + (KI+KPC_re*XL)*Id + (-KPC_im*XL)*Iq
  //    = KPC_im*Vd + KPC_re*Vq + (KI+KPC_re*XL)*Id - KPC_im*XL*Iq
  double KIeff = KI + KPC_re * XL;
  double ve_real = KPC_re * Vd - KPC_im * Vq - KIeff * Iq - KPC_im * XL * Id;
  double ve_imag = KPC_im * Vd + KPC_re * Vq + KIeff * Id  - KPC_im * XL * Iq;
  double VE = sqrt(ve_real * ve_real + ve_imag * ve_imag);
  if (VE < 1e-6) VE = 1e-6;

  // Compute IN and FEX
  double IN  = KC * LadIfd / VE;
  double FEX = computeFex(IN);

  // VB = clamp(VE*FEX, 0, VBMAX)
  VB = std::max(0.0, std::min(VBMAX, VE * FEX));

  // From Efd = VM * VB → VM = Efd / VB
  double VMss = (VB > 1e-6) ? (Efd / VB) : 0.0;

  // Clamp VM to limits
  if (VMss > Vmmax) { Vmmax = VMss; }
  if (VMss < Vmmin) { Vmmin = VMss; }

  // Initialize inner regulator block
  double vrs;
  if (zero_TM) {
    // VMss = KM * vrs → vrs = VMss / KM
    vrs = (fabs(KM) > 1e-6) ? (VMss / KM) : VMss;
    InnerReg_gain_blk.setparams(KM, Vmmin, Vmmax);
  } else {
    InnerReg_blk.setparams(KM, TM, Vmmin, Vmmax, -1000.0, 1000.0);
    vrs = InnerReg_blk.init_given_y(VMss);
  }
  VM = VMss;

  // VG = clamp(KG * Efd, 0, VGMAX)
  double VG = std::max(0.0, std::min(VGMAX, KG * Efd));

  // VR = vrs + VG
  VR = vrs + VG;

  // Widen VR limits if needed
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

  // Initialize regulator block from VR
  double leadlag_in;
  if (zero_TA) {
    leadlag_in = (fabs(KA) > 1e-6) ? (VR / KA) : VR;
  } else {
    leadlag_in = Regulator_blk.init_given_y(VR);
  }

  // Initialize lead-lag
  double vi;
  if (!zero_TB) {
    vi = Leadlag_blk.init_given_y(leadlag_in);
  } else {
    vi = leadlag_in;
  }

  // VI = clamp(Verr, VIMIN, VIMAX) → Verr = vi
  double Verr = vi;
  if (Verr > VIMAX) VIMAX = Verr;
  if (Verr < VIMIN) VIMIN = Verr;

  // Vmeas
  if (!zero_TR) {
    Vmeas = Vmeas_blk.init_given_u(Ec);
  } else {
    Vmeas = Ec;
  }

  Vref = Verr + Vmeas - Vs;
}

void gridpack::dynamic_simulation::Esst3aModel::computeModel(
    double t_inc, IntegrationStage int_flag)
{
  // Voltage measurement
  if (!zero_TR) {
    Vmeas = Vmeas_blk.getoutput(Ec, t_inc, int_flag, true);
  } else {
    Vmeas = Ec;
  }

  // Error + input limiter
  double Verr = Vref + Vs - Vmeas;
  double VI = std::max(VIMIN, std::min(VIMAX, Verr));

  // Lead-lag
  if (!zero_TB) {
    VLL = Leadlag_blk.getoutput(VI, t_inc, int_flag, true);
  } else {
    VLL = VI;
  }

  // Outer voltage regulator
  if (zero_TA) {
    VR = Regulator_gain_blk.getoutput(VLL);
  } else {
    VR = Regulator_blk.getoutput(VLL, t_inc, int_flag, true);
  }

  // Forward-path feedback: VG = clamp(KG * Efd, 0, VGMAX)
  double VG = std::max(0.0, std::min(VGMAX, KG * Efd));
  double vrs = VR - VG;

  // Inner field regulator
  if (zero_TM) {
    VM = InnerReg_gain_blk.getoutput(vrs);
  } else {
    VM = InnerReg_blk.getoutput(vrs, t_inc, int_flag, true);
  }

  // Rectifier terminal voltage (VE)
  double KIeff   = KI + KPC_re * XL;
  double ve_real = KPC_re * Vd - KPC_im * Vq - KIeff * Iq  - KPC_im * XL * Id;
  double ve_imag = KPC_im * Vd + KPC_re * Vq + KIeff * Id  - KPC_im * XL * Iq;
  double VE = sqrt(ve_real * ve_real + ve_imag * ve_imag);
  if (VE < 1e-6) VE = 1e-6;

  // Commutation factor
  double IN  = KC * LadIfd / VE;
  double FEX = computeFex(IN);

  // Rectifier output
  VB = std::max(0.0, std::min(VBMAX, VE * FEX));

  // Output
  Efd = VM * VB;
}

void gridpack::dynamic_simulation::Esst3aModel::predictor(double t_inc, bool flag)
{
  computeModel(t_inc, PREDICTOR);
}

void gridpack::dynamic_simulation::Esst3aModel::corrector(double t_inc, bool flag)
{
  computeModel(t_inc, CORRECTOR);
}

void gridpack::dynamic_simulation::Esst3aModel::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

double gridpack::dynamic_simulation::Esst3aModel::getFieldVoltage()
{
  return Efd;
}

void gridpack::dynamic_simulation::Esst3aModel::setFieldCurrent(double fldc)
{
  LadIfd = fldc;
}

void gridpack::dynamic_simulation::Esst3aModel::setVterminal(double Vm)
{
  Ec = Vm;
  // Default Vq = Vm so VE calculation is non-zero before setVdqIdq is called
  if (Vd == 0.0 && Vq == 0.0) Vq = Vm;
}

void gridpack::dynamic_simulation::Esst3aModel::setVstab(double Vstab)
{
  Vs = Vstab;
}

void gridpack::dynamic_simulation::Esst3aModel::setVdqIdq(
    double vd, double vq, double id, double iq)
{
  Vd = vd;
  Vq = vq;
  Id = id;
  Iq = iq;
}

void gridpack::dynamic_simulation::Esst3aModel::setExtBusNum(int ExtBusNum)
{
  p_bus_num = ExtBusNum;
}

void gridpack::dynamic_simulation::Esst3aModel::setExtGenId(std::string ExtGenId)
{
  p_gen_id = ExtGenId;
}

bool gridpack::dynamic_simulation::Esst3aModel::setState(std::string name,
    double value)
{
  return false;
}

bool gridpack::dynamic_simulation::Esst3aModel::getState(std::string name,
    double *value)
{
  return false;
}
