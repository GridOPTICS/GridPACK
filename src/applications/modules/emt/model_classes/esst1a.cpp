/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   esst1a.cpp
 * @brief  IEEE ST1A static exciter (PSS/E ESST1A) for the EMT module.
 *         Vmeas -> Vi limiter -> two lead-lags -> KA/(1+sTA) -> Va;
 *         Efd = Va - max(0,KLR(Ifd-ILR)) limited to [Vt*VRmin, Vt*VRmax-KC*Ifd].
 */

#include <esst1a.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Esst1a::Esst1a(void)
{
  UEL = VOS = 0;
  TR = VImax = VImin = TC = TB = TC1 = TB1 = 0.0;
  KA = TA = VAmax = VAmin = VRmax = VRmin = 0.0;
  KC = KF = TF = KLR = ILR = 0.0;
  Vmeas = xLL1 = xLL2 = Va = xf = 0.0;
  dVmeas = dxLL1 = dxLL2 = dVa = dxf = 0.0;
  Efd = Ec = Vref = Vs = VF = Ifd = 0.0;
  zero_TR = zero_TB = zero_TB1 = zero_TA = zero_TF = false;
  Vi_at_min = Vi_at_max = Va_at_min = Va_at_max = false;
  nxexc = 5;
}

Esst1a::~Esst1a(void)
{
}

void Esst1a::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxexc = 0;
  *nvar = nxexc;
}

void Esst1a::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTExcModel::load(data,idx);
  if (!data->getValue(EXCITER_UEL, &UEL, idx)) UEL = 0;
  if (!data->getValue(EXCITER_VOS, &VOS, idx)) VOS = 0;
  if (!data->getValue(EXCITER_TR, &TR, idx)) TR = 0.0;
  if (!data->getValue(EXCITER_VIMAX, &VImax, idx)) VImax = 99.0;
  if (!data->getValue(EXCITER_VIMIN, &VImin, idx)) VImin = -99.0;
  if (!data->getValue(EXCITER_TC, &TC, idx)) TC = 0.0;
  if (!data->getValue(EXCITER_TB, &TB, idx)) TB = 0.0;
  if (!data->getValue(EXCITER_TC1, &TC1, idx)) TC1 = 0.0;
  if (!data->getValue(EXCITER_TB1, &TB1, idx)) TB1 = 0.0;
  if (!data->getValue(EXCITER_KA, &KA, idx)) KA = 1.0;
  if (!data->getValue(EXCITER_TA, &TA, idx)) TA = 0.0;
  if (!data->getValue(EXCITER_VAMAX, &VAmax, idx)) VAmax = 99.0;
  if (!data->getValue(EXCITER_VAMIN, &VAmin, idx)) VAmin = -99.0;
  if (!data->getValue(EXCITER_VRMAX, &VRmax, idx)) VRmax = 99.0;
  if (!data->getValue(EXCITER_VRMIN, &VRmin, idx)) VRmin = -99.0;
  if (!data->getValue(EXCITER_KC, &KC, idx)) KC = 0.0;
  if (!data->getValue(EXCITER_KF, &KF, idx)) KF = 0.0;
  if (!data->getValue(EXCITER_TF, &TF, idx)) TF = 0.0;
  if (!data->getValue(EXCITER_KLR, &KLR, idx)) KLR = 0.0;
  if (!data->getValue(EXCITER_ILR, &ILR, idx)) ILR = 0.0;

  zero_TR  = (fabs(TR)  <= 1e-6);
  zero_TB  = (fabs(TB)  <= 1e-6);
  zero_TB1 = (fabs(TB1) <= 1e-6);
  zero_TA  = (fabs(TA)  <= 1e-6);
  zero_TF  = (fabs(TF)  <= 1e-6 || fabs(KF) <= 1e-12);
  if(fabs(KA) < 1e-12) KA = 1.0;

  // Blocks for the explicit path; load() runs before the integration type is set.
  if(!zero_TR) Vmeas_blk.setparams(1.0,TR);
  if(!zero_TB) LL1_blk.setparams(TC,TB);
  if(!zero_TB1) LL2_blk.setparams(TC1,TB1);
  if(!zero_TA) Regulator_blk.setparams(KA,TA,VAmin,VAmax,VAmin,VAmax);
  else Regulator_gain_blk.setparams(KA,VAmin,VAmax);
  if(!zero_TF) {
    double a[2], b[2];
    a[0] = TF; a[1] = 1.0;
    b[0] = KF; b[1] = 0.0;
    Feedback_blk.setcoeffs(a,b);
  }
}

double Esst1a::terminalVoltage()
{
  double vabc[3], vdq0[3];
  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;
  double delta = getGenerator()->getAngle();
  abc2dq0(vabc,p_time,delta,vdq0);
  return sqrt(vdq0[0]*vdq0[0] + vdq0[1]*vdq0[1]);
}

double Esst1a::computeEfd(double Va_in, double Vt, double Ifd_in, bool *clamped)
{
  double e = Va_in - std::max(0.0, KLR*(Ifd_in - ILR));
  double emin = Vt*VRmin, emax = Vt*VRmax - KC*Ifd_in;
  *clamped = (e < emin || e > emax);
  return std::min(std::max(e,emin),emax);
}

void Esst1a::init(gridpack::RealType* xin)
{
  gridpack::RealType *x = xin + offsetb;
  Ec = sqrt(VD*VD + VQ*VQ);
  Efd = getInitialFieldVoltage();
  Ifd = getGenerator()->getFieldCurrent();

  double Vfd = std::max(0.0, KLR*(Ifd - ILR));
  Va = Efd + Vfd;
  double yLL2 = Va/KA, yLL1 = yLL2, Vi = yLL1;
  VF = 0.0;
  Vmeas = Ec;
  Vref = Vi + Vmeas - Vs;

  if(Va < VAmin || Va > VAmax)
    printf("ESST1A bus %d: initial Va=%.4f outside [%.4f, %.4f]\n", busnum, Va, VAmin, VAmax);
  if(Vi < VImin || Vi > VImax)
    printf("ESST1A bus %d: initial Vi=%.4f outside [%.4f, %.4f]\n", busnum, Vi, VImin, VImax);
  if(Efd < Ec*VRmin || Efd > Ec*VRmax - KC*Ifd)
    printf("ESST1A bus %d: initial Efd=%.4f outside [%.4f, %.4f]\n", busnum, Efd, Ec*VRmin, Ec*VRmax - KC*Ifd);

  if(integrationtype != IMPLICIT) {
    if(!zero_TR) Vmeas_blk.init_given_u(Ec);
    if(!zero_TF) Feedback_blk.init_given_u(Efd);
    if(!zero_TB) LL1_blk.init_given_u(Vi);
    if(!zero_TB1) LL2_blk.init_given_u(yLL1);
    if(!zero_TA) Regulator_blk.init_given_y(Va);
  } else {
    xf   = zero_TF  ? 0.0 : -KF/TF*Efd;
    xLL1 = zero_TB  ? Vi   : (1.0 - TC/TB)*Vi;
    xLL2 = zero_TB1 ? yLL1 : (1.0 - TC1/TB1)*yLL1;
    x[0] = Vmeas; x[1] = xLL1; x[2] = xLL2; x[3] = Va; x[4] = xf;
  }
}

void Esst1a::preStep(double time, double timestep)
{
  if(integrationtype != EXPLICIT) return;
  Ec = terminalVoltage();
  Vmeas = zero_TR ? Ec : Vmeas_blk.getoutput(Ec,timestep,true);
  VF = zero_TF ? 0.0 : Feedback_blk.getoutput(Efd,timestep,true);
  double Vi = std::min(std::max(Vref - Vmeas + Vs - VF, VImin), VImax);
  double yLL1 = zero_TB  ? Vi   : LL1_blk.getoutput(Vi,timestep,true);
  double yLL2 = zero_TB1 ? yLL1 : LL2_blk.getoutput(yLL1,timestep,true);
  Va = zero_TA ? Regulator_gain_blk.getoutput(yLL2) : Regulator_blk.getoutput(yLL2,timestep,true);
  Ifd = getGenerator()->getFieldCurrent();
  bool clamped;
  Efd = computeEfd(Va, Ec, Ifd, &clamped);
}

void Esst1a::postStep(double time)
{
}

bool Esst1a::serialWrite(char *string, const int bufsize, const char *signal)
{
  return false;
}

void Esst1a::write(const char* signal, char* string)
{
}

void Esst1a::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val + offsetb;
  if(integrationtype == EXPLICIT) return;
  if(p_mode == XVECTOBUS) {
    Vmeas = values[0]; xLL1 = values[1]; xLL2 = values[2]; Va = values[3]; xf = values[4];
  } else if(p_mode == XDOTVECTOBUS) {
    dVmeas = values[0]; dxLL1 = values[1]; dxLL2 = values[2]; dVa = values[3]; dxf = values[4];
  }
}

void Esst1a::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values + offsetb;
  if(integrationtype == EXPLICIT || p_mode != RESIDUAL_EVAL) return;

  Ec = terminalVoltage();
  Ifd = getGenerator()->getFieldCurrent();
  bool clamped;
  Efd = computeEfd(Va, Ec, Ifd, &clamped);

  f[0] = zero_TR ? (-Vmeas + Ec) : (-Vmeas + Ec)/TR - dVmeas;

  double Vf = zero_TF ? 0.0 : xf + KF/TF*Efd;
  double Vi = Vref - Vmeas + Vs - Vf;
  if(Vi_at_max) Vi = VImax;
  else if(Vi_at_min) Vi = VImin;

  double yLL1, yLL2;
  if(zero_TB) { f[1] = -xLL1 + Vi; yLL1 = xLL1; }
  else { f[1] = (-xLL1 + (1.0 - TC/TB)*Vi)/TB - dxLL1; yLL1 = xLL1 + TC/TB*Vi; }

  if(zero_TB1) { f[2] = -xLL2 + yLL1; yLL2 = xLL2; }
  else { f[2] = (-xLL2 + (1.0 - TC1/TB1)*yLL1)/TB1 - dxLL2; yLL2 = xLL2 + TC1/TB1*yLL1; }

  if(Va_at_min) f[3] = Va - VAmin;
  else if(Va_at_max) f[3] = Va - VAmax;
  else f[3] = zero_TA ? (-Va + KA*yLL2) : (-Va + KA*yLL2)/TA - dVa;

  f[4] = zero_TF ? -xf : (-xf - KF/TF*Efd)/TF - dxf;
}

int Esst1a::matrixNumValues()
{
  return (integrationtype == IMPLICIT) ? 21 : 0;
}

void Esst1a::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  if(integrationtype != IMPLICIT) { *nvals = 0; return; }

  int Vmeas_idx = p_gloc, xLL1_idx = p_gloc+1, xLL2_idx = p_gloc+2, Va_idx = p_gloc+3, xf_idx = p_gloc+4;
  int delta_idx;
  double delta = getGenerator()->getAngle(&delta_idx);

  // Ec = |V_dq| partials w.r.t. va,vb,vc and delta
  double Tdq0[3][3], dTdq0ddelta[3][3], vabc[3], vdq0[3];
  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;
  getTdq0(p_time,delta,Tdq0);
  getdTdq0dtheta(p_time,delta,dTdq0ddelta);
  abc2dq0(vabc,p_time,delta,vdq0);
  double Vd = vdq0[0], Vq = vdq0[1];
  double Ecv = sqrt(Vd*Vd + Vq*Vq);
  double dEc_dVd = (Ecv > 1e-12) ? Vd/Ecv : 0.0, dEc_dVq = (Ecv > 1e-12) ? Vq/Ecv : 0.0;
  double dEc_dv[3], dEc_ddelta;
  for(int j = 0; j < 3; j++) dEc_dv[j] = dEc_dVd*Tdq0[0][j] + dEc_dVq*Tdq0[1][j];
  dEc_ddelta = dEc_dVd*(dTdq0ddelta[0][0]*vabc[0] + dTdq0ddelta[0][1]*vabc[1] + dTdq0ddelta[0][2]*vabc[2])
             + dEc_dVq*(dTdq0ddelta[1][0]*vabc[0] + dTdq0ddelta[1][1]*vabc[1] + dTdq0ddelta[1][2]*vabc[2]);

  // Row Vmeas
  double s0 = zero_TR ? 1.0 : 1.0/TR;
  rows[ctr] = Vmeas_idx; cols[ctr] = Vmeas_idx; values[ctr++] = -s0 - (zero_TR ? 0.0 : shift);
  rows[ctr] = Vmeas_idx; cols[ctr] = delta_idx; values[ctr++] = s0*dEc_ddelta;
  for(int j = 0; j < 3; j++) { rows[ctr] = Vmeas_idx; cols[ctr] = p_glocvoltage+j; values[ctr++] = s0*dEc_dv[j]; }

  // Vi = Vref - Vmeas + Vs - Vf, Vf = xf + KF/TF*Efd, Efd = Va (unless clamped)
  bool clamped;
  computeEfd(Va, Ecv, Ifd, &clamped);
  double dEfd_dVa = clamped ? 0.0 : 1.0;
  double dVi_dVmeas = -1.0, dVi_dxf = zero_TF ? 0.0 : -1.0, dVi_dVa = zero_TF ? 0.0 : -KF/TF*dEfd_dVa;
  if(Vi_at_min || Vi_at_max) dVi_dVmeas = dVi_dxf = dVi_dVa = 0.0;

  // Row xLL1; yLL1 partials
  double c1 = zero_TB ? 1.0 : (1.0 - TC/TB)/TB;
  rows[ctr] = xLL1_idx; cols[ctr] = xLL1_idx;  values[ctr++] = zero_TB ? -1.0 : -1.0/TB - shift;
  rows[ctr] = xLL1_idx; cols[ctr] = Vmeas_idx; values[ctr++] = c1*dVi_dVmeas;
  rows[ctr] = xLL1_idx; cols[ctr] = xf_idx;    values[ctr++] = c1*dVi_dxf;
  rows[ctr] = xLL1_idx; cols[ctr] = Va_idx;    values[ctr++] = c1*dVi_dVa;
  double dy1_dxLL1 = 1.0, dy1_dVi = zero_TB ? 0.0 : TC/TB;

  // Row xLL2; yLL2 partials
  double c2 = zero_TB1 ? 1.0 : (1.0 - TC1/TB1)/TB1;
  rows[ctr] = xLL2_idx; cols[ctr] = xLL2_idx;  values[ctr++] = zero_TB1 ? -1.0 : -1.0/TB1 - shift;
  rows[ctr] = xLL2_idx; cols[ctr] = xLL1_idx;  values[ctr++] = c2*dy1_dxLL1;
  rows[ctr] = xLL2_idx; cols[ctr] = Vmeas_idx; values[ctr++] = c2*dy1_dVi*dVi_dVmeas;
  rows[ctr] = xLL2_idx; cols[ctr] = xf_idx;    values[ctr++] = c2*dy1_dVi*dVi_dxf;
  rows[ctr] = xLL2_idx; cols[ctr] = Va_idx;    values[ctr++] = c2*dy1_dVi*dVi_dVa;
  double dy2_dxLL2 = 1.0, dy2_dy1 = zero_TB1 ? 0.0 : TC1/TB1;

  // Row Va
  if(Va_at_min || Va_at_max) {
    rows[ctr] = Va_idx; cols[ctr] = Va_idx;    values[ctr++] = 1.0;
    rows[ctr] = Va_idx; cols[ctr] = xLL2_idx;  values[ctr++] = 0.0;
    rows[ctr] = Va_idx; cols[ctr] = xLL1_idx;  values[ctr++] = 0.0;
    rows[ctr] = Va_idx; cols[ctr] = Vmeas_idx; values[ctr++] = 0.0;
    rows[ctr] = Va_idx; cols[ctr] = xf_idx;    values[ctr++] = 0.0;
  } else {
    double k = zero_TA ? KA : KA/TA;
    double dchain = dy2_dy1*dy1_dVi;   // d(yLL2)/d(Vi)
    rows[ctr] = Va_idx; cols[ctr] = Va_idx;    values[ctr++] = (zero_TA ? -1.0 : -1.0/TA - shift) + k*dchain*dVi_dVa;
    rows[ctr] = Va_idx; cols[ctr] = xLL2_idx;  values[ctr++] = k*dy2_dxLL2;
    rows[ctr] = Va_idx; cols[ctr] = xLL1_idx;  values[ctr++] = k*dy2_dy1*dy1_dxLL1;
    rows[ctr] = Va_idx; cols[ctr] = Vmeas_idx; values[ctr++] = k*dchain*dVi_dVmeas;
    rows[ctr] = Va_idx; cols[ctr] = xf_idx;    values[ctr++] = k*dchain*dVi_dxf;
  }

  // Row xf
  rows[ctr] = xf_idx; cols[ctr] = xf_idx; values[ctr++] = zero_TF ? -1.0 : -1.0/TF - shift;
  rows[ctr] = xf_idx; cols[ctr] = Va_idx; values[ctr++] = zero_TF ? 0.0 : -(KF/TF)/TF*dEfd_dVa;

  *nvals = ctr;
}

double Esst1a::getFieldVoltage()
{
  if(integrationtype == IMPLICIT) {
    bool clamped;
    Efd = computeEfd(Va, terminalVoltage(), Ifd, &clamped);
  }
  return Efd;
}

double Esst1a::getFieldVoltage(int *Efd_gloc)
{
  *Efd_gloc = (integrationtype == IMPLICIT) ? p_gloc + 3 : -1;
  return getFieldVoltage();
}

bool Esst1a::getFieldVoltagePartialDerivatives(int *xexc_loc, double *dEfd_dxexc, double *dEfd_dxgen)
{
  return false;
}

void Esst1a::resetEventFlags()
{
  Vi_at_min = Vi_at_max = Va_at_min = Va_at_max = false;
}

void Esst1a::eventFunction(const double&t, gridpack::RealType *state, std::vector<gridpack::RealType >& evalues)
{
  int off = getLocalOffset();
  Vmeas = state[off]; xLL1 = state[off+1]; xLL2 = state[off+2]; Va = state[off+3]; xf = state[off+4];

  double Vf = zero_TF ? 0.0 : xf + KF/TF*Efd;
  double Vi = Vref - Vmeas + Vs - Vf;
  double yLL1 = zero_TB  ? xLL1 : xLL1 + TC/TB*Vi;
  double yLL2 = zero_TB1 ? xLL2 : xLL2 + TC1/TB1*yLL1;
  double dVa_dt = zero_TA ? 0.0 : (-Va + KA*yLL2)/TA;

  evalues[0] = Vi_at_min ? (VImin - Vi) : (Vi - VImin);
  evalues[1] = Vi_at_max ? (Vi - VImax) : (VImax - Vi);
  evalues[2] = Va_at_min ? -dVa_dt : (Va - VAmin);
  evalues[3] = Va_at_max ?  dVa_dt : (VAmax - Va);
}

void Esst1a::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int off = getLocalOffset();
  Vmeas = state[off]; xLL1 = state[off+1]; xLL2 = state[off+2]; Va = state[off+3]; xf = state[off+4];

  double Vf = zero_TF ? 0.0 : xf + KF/TF*Efd;
  double Vi = Vref - Vmeas + Vs - Vf;
  double yLL1 = zero_TB  ? xLL1 : xLL1 + TC/TB*Vi;
  double yLL2 = zero_TB1 ? xLL2 : xLL2 + TC1/TB1*yLL1;
  double dVa_dt = zero_TA ? 0.0 : (-Va + KA*yLL2)/TA;

  if(triggered[0]) Vi_at_min = !Vi_at_min;
  if(triggered[1]) Vi_at_max = !Vi_at_max;
  if(triggered[2]) Va_at_min = (!Va_at_min && dVa_dt < 0);
  if(triggered[3]) Va_at_max = (!Va_at_max && dVa_dt > 0);
}

void Esst1a::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  if(integrationtype == IMPLICIT) {
    gridpack::math::RealDAESolver::EventPtr e(new Esst1aEvent(this));
    eman->add(e);
  }
}

void Esst1aEvent::p_update(const double& t, gridpack::RealType *state)
{
  p_exc->eventFunction(t,state,p_current);
}

void Esst1aEvent::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_exc->eventHandlerFunction(triggered,t,state);
}
