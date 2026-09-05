/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   gensal.cpp
 * @brief  Salient-pole machine (PSS/E GENSAL) for the EMT module. States:
 *         psid psiq psi0 | Eqp psi1d psi2q delta dw | ia ib ic.
 */

#include <gensal.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Gensal::Gensal(void)
{
  nxgen = 11;
  flux_speed_sensitivity = 1;
  enableSat = false;
  TM = 0.0;
  LadIfd = 0.0;
}

Gensal::~Gensal(void)
{
}

void Gensal::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgen = 6;
  *nvar = nxgen;
}

void Gensal::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGenModel::load(data,idx);

  gridpack::ComplexType Zsource;
  data->getValue(BUS_NUMBER, &bid);
  if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H, &H, idx)) H = 0.0;
  if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0, &D, idx)) D = 0.0;
  data->getValue(GENERATOR_ZSOURCE,&Zsource,idx);
  Ra = real(Zsource);
  if (!data->getValue(GENERATOR_RESISTANCE, &Ra, idx)) Ra = 0.0;
  if (!data->getValue(GENERATOR_XD, &Xd, idx)) Xd = 0.0;
  if (!data->getValue(GENERATOR_XQ, &Xq, idx)) Xq = 0.0;
  if (!data->getValue(GENERATOR_XDP, &Xdp, idx)) Xdp = 0.0;
  if (!data->getValue(GENERATOR_XDPP, &Xdpp, idx)) Xdpp = 0.0;
  if (!data->getValue(GENERATOR_XL, &Xl, idx)) Xl = 0.0;
  if (!data->getValue(GENERATOR_TDOP, &Tdop, idx)) Tdop = 0.0;
  if (!data->getValue(GENERATOR_TDOPP, &Tdopp, idx)) Tdopp = 0.0;
  if (!data->getValue(GENERATOR_TQOPP, &Tqopp, idx)) Tqopp = 0.0;
  if (!data->getValue(GENERATOR_S1, &S10, idx)) S10 = 0.0;
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12 = 0.0;

  enableSat = (fabs(S10*S12) >= 1e-6);

  // Parameter guards (same rules as Genrou / PowerWorld).
  if (Xdp > Xd) {
    printf("GENSAL bus %d: Auto-correct Xd'=%.6f > Xd=%.6f -> Xd'=%.6f\n",
           bid, Xdp, Xd, 0.8*Xd);
    Xdp = 0.8*Xd;
  }
  if (Xdpp > Xdp) {
    printf("GENSAL bus %d: Auto-correct Xd''=%.6f > Xd'=%.6f -> Xd''=%.6f\n",
           bid, Xdpp, Xdp, 0.8*Xdp);
    Xdpp = 0.8*Xdp;
  }
  if (Xl >= Xdpp) {
    printf("GENSAL bus %d: Auto-correct Xl=%.6f >= Xd''=%.6f -> Xl=%.6f\n",
           bid, Xl, Xdpp, 0.8*Xdpp);
    Xl = 0.8*Xdpp;
  }
  if (Xq < Xdpp) {
    printf("GENSAL bus %d: Auto-correct Xq=%.6f < Xd''=%.6f -> Xq=%.6f\n",
           bid, Xq, Xdpp, Xdpp);
    Xq = Xdpp;
  }
  if (Tdop <= 0.0) Tdop = 1e-3;
  if (Tdopp <= 0.0) Tdopp = 1e-3;
  if (Tqopp <= 0.0) Tqopp = 1e-3;
}

// Scaled quadratic saturation Se(x) = B(x-A)^2/x fitted to S(1.0), S(1.2).
double Gensal::Sat(double x)
{
  if (!enableSat || x <= 1e-6) return 0.0;
  double R = 1.2*S12/S10;
  double sqrtR = sqrt(R);
  double A = (1.2 - sqrtR)/(1.0 - sqrtR);
  double B = S10/((1.0 - A)*(1.0 - A));
  double tmp = x - A;
  if (tmp < 0.0) return 0.0;
  return B*tmp*tmp/x;
}

double Gensal::dSat(double x)
{
  if (!enableSat || x <= 1e-6) return 0.0;
  double R = 1.2*S12/S10;
  double sqrtR = sqrt(R);
  double A = (1.2 - sqrtR)/(1.0 - sqrtR);
  double B = S10/((1.0 - A)*(1.0 - A));
  if (x - A < 0.0) return 0.0;
  return B*(x - A)*(x + A)/(x*x);
}

void Gensal::init(gridpack::RealType* xin)
{
  gridpack::RealType *x = xin + offsetb;
  double Pg = pg/mbase, Qg = qg/mbase;

  VD = p_Vm0*cos(p_Va0);
  VQ = p_Vm0*sin(p_Va0);
  gridpack::ComplexType V(VD,VQ), S(Pg,Qg);
  gridpack::ComplexType I = conj(S/V);
  double Im = abs(I), Ia = arg(I);

  iabc[0] = Im*sin(OMEGA_S*p_time + Ia);
  iabc[1] = Im*sin(OMEGA_S*p_time + Ia - TWOPI_OVER_THREE);
  iabc[2] = Im*sin(OMEGA_S*p_time + Ia + TWOPI_OVER_THREE);
  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;

  // No q-axis saturation in GENSAL: E = V + (Ra + jXq) I lies on the q axis.
  delta = arg(V + gridpack::ComplexType(Ra,Xq)*I);
  double theta = delta - PI/2.0;
  abc2dq0(vabc,p_time,theta,vdq0);
  abc2dq0(iabc,p_time,theta,idq0);
  double Vd = vdq0[0], Vq = vdq0[1];
  double Id = idq0[0], Iq = idq0[1];

  psid = Ra*Iq + Vq;
  psiq = -Ra*Id - Vd;
  psi0 = 0.0;
  dw = 0.0;

  psi1d = psid + Xl*Id;
  Eqp   = psi1d + (Xdp - Xl)*Id;
  psi2q = psiq + Xdpp*Iq;          // q-axis subtransient flux

  TM = psid*Iq - psiq*Id;
  LadIfd = Eqp*(1.0 + Sat(Eqp)) + (Xd - Xdp)*Id;
  Efd = LadIfd;

  double scal = mbase/sbase;
  if(integrationtype != EXPLICIT) {
    x[0] = psid; x[1] = psiq; x[2] = psi0;
    x[3] = Eqp;  x[4] = psi1d; x[5] = psi2q;
    x[6] = delta; x[7] = dw;
    x[8] = iabc[0]*scal; x[9] = iabc[1]*scal; x[10] = iabc[2]*scal;
  } else {
    x[0] = psid; x[1] = psiq; x[2] = psi0;
    x[3] = iabc[0]*scal; x[4] = iabc[1]*scal; x[5] = iabc[2]*scal;
  }
}

bool Gensal::serialWrite(char *string, const int bufsize,const char *signal)
{
  if(!strcmp(signal,"header")) {
    sprintf(string,", %d_%s_V,%d_%s_Pg,%d_%s_Qg,%d_%s_delta, %d_%s_dw",
            busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str(),
            busnum,id.c_str(),busnum,id.c_str());
    return true;
  } else if(!strcmp(signal,"monitor")) {
    double Vm, Pgen, Qgen, dspd;
    if(p_online) {
      Vm = sqrt(vdq0[0]*vdq0[0] + vdq0[1]*vdq0[1]);
      dspd = dw;
    } else {
      Vm = p_Vm0;
      dspd = 0.0;
    }
    Pgen = p_online*(vdq0[0]*idq0[0] + vdq0[1]*idq0[1])*mbase/sbase;
    Qgen = p_online*(vdq0[1]*idq0[0] - vdq0[0]*idq0[1])*mbase/sbase;
    sprintf(string,", %.17g,%.17g,%.17g,%.17g,%.17g",Vm,Pgen,Qgen,delta,dspd);
    return true;
  }
  return false;
}

void Gensal::write(const char* signal, char* string)
{
}

void Gensal::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values + offsetb;
  if(integrationtype == EXPLICIT) {
    if(p_mode == XVECTOBUS) {
      psid = x[0]; psiq = x[1]; psi0 = x[2];
      iabc[0] = x[3]; iabc[1] = x[4]; iabc[2] = x[5];
    } else if(p_mode == XDOTVECTOBUS) {
      dpsid = x[0]; dpsiq = x[1]; dpsi0 = x[2];
      diabc[0] = x[3]; diabc[1] = x[4]; diabc[2] = x[5];
    }
  } else {
    if(p_mode == XVECTOBUS) {
      psid = x[0]; psiq = x[1]; psi0 = x[2];
      Eqp = x[3]; psi1d = x[4]; psi2q = x[5];
      delta = x[6]; dw = x[7];
      iabc[0] = x[8]; iabc[1] = x[9]; iabc[2] = x[10];
    } else if(p_mode == XDOTVECTOBUS) {
      dpsid = x[0]; dpsiq = x[1]; dpsi0 = x[2];
      dEqp = x[3]; dpsi1d = x[4]; dpsi2q = x[5];
      ddelta = x[6]; ddw = x[7];
      diabc[0] = x[8]; diabc[1] = x[9]; diabc[2] = x[10];
    }
  }
}

void Gensal::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values + offsetb;
  if(p_mode != RESIDUAL_EVAL) return;

  double theta = delta - PI/2.0;
  double tempd1 = (Xdpp - Xl)/(Xdp - Xl);
  double tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
  double param1 = (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));

  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;
  abc2dq0(vabc,p_time,theta,vdq0);
  double Vd = vdq0[0], Vq = vdq0[1], V0 = vdq0[2];

  double Id = (psid - tempd1*Eqp - tempd2*psi1d)/-Xdpp;
  double Iq = (psiq - psi2q)/-Xdpp;
  double I0 = psi0/-Xl;
  idq0[0] = Id; idq0[1] = Iq; idq0[2] = I0;

  if(hasExciter()) Efd = getExciter()->getFieldVoltage();
  if(hasGovernor()) TM = getGovernor()->getMechanicalPower();

  double dpsi1ddt = -psi1d + Eqp - (Xdp - Xl)*Id;
  LadIfd = Eqp*(1.0 + Sat(Eqp)) + (Xd - Xdp)*(Id + param1*dpsi1ddt);

  double igen[3];
  dq02abc(idq0,p_time,theta,igen);
  double scal = mbase/sbase;
  double k = flux_speed_sensitivity;

  if(integrationtype != EXPLICIT) {
    f[0] = OMEGA_S*(Ra*Id + (1 + k*dw)*psiq + Vd) - dpsid;
    f[1] = OMEGA_S*(Ra*Iq - (1 + k*dw)*psid + Vq) - dpsiq;
    f[2] = OMEGA_S*(Ra*I0 + V0) - dpsi0;
    f[3] = (Efd - LadIfd)/Tdop - dEqp;
    f[4] = dpsi1ddt/Tdopp - dpsi1d;
    f[5] = (-psi2q - (Xq - Xdpp)*Iq)/Tqopp - dpsi2q;
    f[6] = dw*OMEGA_S - ddelta;
    f[7] = 1.0/(2.0*H)*((TM - D*dw)/(1 + dw) - (psid*Iq - psiq*Id)) - ddw;
    f[8]  = igen[0]*scal - iabc[0];
    f[9]  = igen[1]*scal - iabc[1];
    f[10] = igen[2]*scal - iabc[2];
  } else if(p_online) {
    f[0] = OMEGA_S*(Ra*Id + (1 + k*dw)*psiq + Vd) - dpsid;
    f[1] = OMEGA_S*(Ra*Iq - (1 + k*dw)*psid + Vq) - dpsiq;
    f[2] = OMEGA_S*(Ra*I0 + V0) - dpsi0;
    f[3] = igen[0]*scal - iabc[0];
    f[4] = igen[1]*scal - iabc[1];
    f[5] = igen[2]*scal - iabc[2];
  } else {
    f[0] = psid; f[1] = psiq; f[2] = psi0;
    f[3] = iabc[0]; f[4] = iabc[1]; f[5] = iabc[2];
  }
}

void Gensal::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = p_online*iabc[0];
  *ib = p_online*iabc[1];
  *ic = p_online*iabc[2];
}

void Gensal::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = (integrationtype != EXPLICIT) ? p_gloc + 8 : p_gloc + 3;
}

int Gensal::matrixNumValues()
{
  int numVals = 0;
  if(integrationtype != EXPLICIT) {
    numVals = 62;
    if(hasExciter()) numVals += 1;
    if(hasGovernor()) numVals += 1;
  } else {
    numVals = 27;
  }
  return numVals;
}

void Gensal::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0, j;
  double k = flux_speed_sensitivity;
  double Tdq0[3][3], Tdq0inv[3][3];
  double theta = delta - PI/2.0;
  double scal = mbase/sbase;
  double dId_dpsid = -1/Xdpp, dIq_dpsiq = -1/Xdpp, dI0_dpsi0 = -1/Xl;

  getTdq0(p_time,theta,Tdq0);
  getTdq0inv(p_time,theta,Tdq0inv);

  if(integrationtype == EXPLICIT) {
    int psid_idx = p_gloc, psiq_idx = p_gloc+1, psi0_idx = p_gloc+2;
    int i_idx[3] = { p_gloc+3, p_gloc+4, p_gloc+5 };
    if(!p_online) {
      for(j = 0; j < 6; j++) { rows[ctr] = cols[ctr] = p_gloc + j; values[ctr] = 1.0; ctr++; }
      *nvals = ctr;
      return;
    }
    rows[ctr] = psid_idx; cols[ctr] = psid_idx; values[ctr++] = OMEGA_S*Ra*dId_dpsid - shift;
    rows[ctr] = psid_idx; cols[ctr] = psiq_idx; values[ctr++] = OMEGA_S*(1 + k*dw);
    for(j = 0; j < 3; j++) { rows[ctr] = psid_idx; cols[ctr] = p_glocvoltage+j; values[ctr++] = OMEGA_S*Tdq0[0][j]; }
    rows[ctr] = psiq_idx; cols[ctr] = psid_idx; values[ctr++] = -OMEGA_S*(1 + k*dw);
    rows[ctr] = psiq_idx; cols[ctr] = psiq_idx; values[ctr++] = OMEGA_S*Ra*dIq_dpsiq - shift;
    for(j = 0; j < 3; j++) { rows[ctr] = psiq_idx; cols[ctr] = p_glocvoltage+j; values[ctr++] = OMEGA_S*Tdq0[1][j]; }
    rows[ctr] = psi0_idx; cols[ctr] = psi0_idx; values[ctr++] = OMEGA_S*Ra*dI0_dpsi0 - shift;
    for(j = 0; j < 3; j++) { rows[ctr] = psi0_idx; cols[ctr] = p_glocvoltage+j; values[ctr++] = OMEGA_S*Tdq0[2][j]; }
    for(int r = 0; r < 3; r++) {
      rows[ctr] = i_idx[r]; cols[ctr] = psid_idx; values[ctr++] = scal*Tdq0inv[r][0]*dId_dpsid;
      rows[ctr] = i_idx[r]; cols[ctr] = psiq_idx; values[ctr++] = scal*Tdq0inv[r][1]*dIq_dpsiq;
      rows[ctr] = i_idx[r]; cols[ctr] = psi0_idx; values[ctr++] = scal*Tdq0inv[r][2]*dI0_dpsi0;
      rows[ctr] = i_idx[r]; cols[ctr] = i_idx[r]; values[ctr++] = -1.0;
    }
    *nvals = ctr;
    return;
  }

  int psid_idx = p_gloc, psiq_idx = p_gloc+1, psi0_idx = p_gloc+2;
  int Eqp_idx = p_gloc+3, psi1d_idx = p_gloc+4, psi2q_idx = p_gloc+5;
  int delta_idx = p_gloc+6, dw_idx = p_gloc+7;
  int i_idx[3] = { p_gloc+8, p_gloc+9, p_gloc+10 };

  double tempd1 = (Xdpp - Xl)/(Xdp - Xl);
  double tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
  double param1 = (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));
  double dId_dEqp = tempd1/Xdpp, dId_dpsi1d = tempd2/Xdpp;
  double dIq_dpsi2q = 1/Xdpp;

  double dTdq0_ddelta[3][3], dTdq0inv_ddelta[3][3], dvdq0_ddelta[3];
  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;
  getdTdq0dtheta(p_time,theta,dTdq0_ddelta);
  getdTdq0invdtheta(p_time,theta,dTdq0inv_ddelta);
  matvecmult3x3(dTdq0_ddelta,vabc,dvdq0_ddelta);

  double Id = (psid - tempd1*Eqp - tempd2*psi1d)/-Xdpp;
  double Iq = (psiq - psi2q)/-Xdpp;
  double I0 = psi0/-Xl;

  // Stator flux rows
  rows[ctr] = psid_idx; cols[ctr] = psid_idx;  values[ctr++] = OMEGA_S*Ra*dId_dpsid - shift;
  rows[ctr] = psid_idx; cols[ctr] = psiq_idx;  values[ctr++] = OMEGA_S*(1 + k*dw);
  rows[ctr] = psid_idx; cols[ctr] = Eqp_idx;   values[ctr++] = OMEGA_S*Ra*dId_dEqp;
  rows[ctr] = psid_idx; cols[ctr] = psi1d_idx; values[ctr++] = OMEGA_S*Ra*dId_dpsi1d;
  rows[ctr] = psid_idx; cols[ctr] = delta_idx; values[ctr++] = OMEGA_S*dvdq0_ddelta[0];
  rows[ctr] = psid_idx; cols[ctr] = dw_idx;    values[ctr++] = OMEGA_S*k*psiq;
  for(j = 0; j < 3; j++) { rows[ctr] = psid_idx; cols[ctr] = p_glocvoltage+j; values[ctr++] = OMEGA_S*Tdq0[0][j]; }

  rows[ctr] = psiq_idx; cols[ctr] = psid_idx;  values[ctr++] = -OMEGA_S*(1 + k*dw);
  rows[ctr] = psiq_idx; cols[ctr] = psiq_idx;  values[ctr++] = OMEGA_S*Ra*dIq_dpsiq - shift;
  rows[ctr] = psiq_idx; cols[ctr] = psi2q_idx; values[ctr++] = OMEGA_S*Ra*dIq_dpsi2q;
  rows[ctr] = psiq_idx; cols[ctr] = delta_idx; values[ctr++] = OMEGA_S*dvdq0_ddelta[1];
  rows[ctr] = psiq_idx; cols[ctr] = dw_idx;    values[ctr++] = -OMEGA_S*k*psid;
  for(j = 0; j < 3; j++) { rows[ctr] = psiq_idx; cols[ctr] = p_glocvoltage+j; values[ctr++] = OMEGA_S*Tdq0[1][j]; }

  rows[ctr] = psi0_idx; cols[ctr] = psi0_idx;  values[ctr++] = OMEGA_S*Ra*dI0_dpsi0 - shift;
  rows[ctr] = psi0_idx; cols[ctr] = delta_idx; values[ctr++] = OMEGA_S*dvdq0_ddelta[2];
  for(j = 0; j < 3; j++) { rows[ctr] = psi0_idx; cols[ctr] = p_glocvoltage+j; values[ctr++] = OMEGA_S*Tdq0[2][j]; }

  // Eqp row: f = (Efd - LadIfd)/Tdop, LadIfd = Eqp(1+Se) + (Xd-Xdp)(Id + param1*dpsi1ddt)
  double dpsi1ddt_dpsid  = -(Xdp - Xl)*dId_dpsid;
  double dpsi1ddt_dEqp   = 1.0 - (Xdp - Xl)*dId_dEqp;
  double dpsi1ddt_dpsi1d = -1.0 - (Xdp - Xl)*dId_dpsi1d;
  double dLad_dpsid  = (Xd - Xdp)*(dId_dpsid + param1*dpsi1ddt_dpsid);
  double dLad_dEqp   = 1.0 + Sat(Eqp) + Eqp*dSat(Eqp) + (Xd - Xdp)*(dId_dEqp + param1*dpsi1ddt_dEqp);
  double dLad_dpsi1d = (Xd - Xdp)*(dId_dpsi1d + param1*dpsi1ddt_dpsi1d);
  rows[ctr] = Eqp_idx; cols[ctr] = psid_idx;  values[ctr++] = -dLad_dpsid/Tdop;
  rows[ctr] = Eqp_idx; cols[ctr] = Eqp_idx;   values[ctr++] = -dLad_dEqp/Tdop - shift;
  rows[ctr] = Eqp_idx; cols[ctr] = psi1d_idx; values[ctr++] = -dLad_dpsi1d/Tdop;
  if(hasExciter()) {
    int Efd_idx;
    getExciter()->getFieldVoltage(&Efd_idx);
    if(Efd_idx >= 0) { rows[ctr] = Eqp_idx; cols[ctr] = Efd_idx; values[ctr++] = 1.0/Tdop; }
  }

  // psi1d row
  rows[ctr] = psi1d_idx; cols[ctr] = psid_idx;  values[ctr++] = dpsi1ddt_dpsid/Tdopp;
  rows[ctr] = psi1d_idx; cols[ctr] = Eqp_idx;   values[ctr++] = dpsi1ddt_dEqp/Tdopp;
  rows[ctr] = psi1d_idx; cols[ctr] = psi1d_idx; values[ctr++] = dpsi1ddt_dpsi1d/Tdopp - shift;

  // psi2q row: f = (-psi2q - (Xq-Xdpp)*Iq)/Tqopp
  rows[ctr] = psi2q_idx; cols[ctr] = psiq_idx;  values[ctr++] = -(Xq - Xdpp)*dIq_dpsiq/Tqopp;
  rows[ctr] = psi2q_idx; cols[ctr] = psi2q_idx; values[ctr++] = (-1.0 - (Xq - Xdpp)*dIq_dpsi2q)/Tqopp - shift;

  // delta row
  rows[ctr] = delta_idx; cols[ctr] = delta_idx; values[ctr++] = -shift;
  rows[ctr] = delta_idx; cols[ctr] = dw_idx;    values[ctr++] = OMEGA_S;

  // dw row: f = Minv*((TM - D*dw)/(1+dw) - Te), Te = psid*Iq - psiq*Id
  double Minv = 1.0/(2.0*H);
  int TM_idx = -1;
  if(hasGovernor()) TM = getGovernor()->getMechanicalPower(&TM_idx);
  rows[ctr] = dw_idx; cols[ctr] = psid_idx;  values[ctr++] = -Minv*(Iq - psiq*dId_dpsid);
  rows[ctr] = dw_idx; cols[ctr] = psiq_idx;  values[ctr++] = -Minv*(psid*dIq_dpsiq - Id);
  rows[ctr] = dw_idx; cols[ctr] = Eqp_idx;   values[ctr++] = -Minv*(-psiq*dId_dEqp);
  rows[ctr] = dw_idx; cols[ctr] = psi1d_idx; values[ctr++] = -Minv*(-psiq*dId_dpsi1d);
  rows[ctr] = dw_idx; cols[ctr] = psi2q_idx; values[ctr++] = -Minv*(psid*dIq_dpsi2q);
  rows[ctr] = dw_idx; cols[ctr] = dw_idx;
  values[ctr++] = Minv*(-D/(1 + dw) - (TM - D*dw)/((1 + dw)*(1 + dw))) - shift;
  if(hasGovernor() && TM_idx >= 0) {
    rows[ctr] = dw_idx; cols[ctr] = TM_idx; values[ctr++] = Minv/(1 + dw);
  }

  // Current rows: igen = Tdq0inv * idq0 (machine base) scaled to system base
  for(int r = 0; r < 3; r++) {
    rows[ctr] = i_idx[r]; cols[ctr] = psid_idx;  values[ctr++] = scal*Tdq0inv[r][0]*dId_dpsid;
    rows[ctr] = i_idx[r]; cols[ctr] = psiq_idx;  values[ctr++] = scal*Tdq0inv[r][1]*dIq_dpsiq;
    rows[ctr] = i_idx[r]; cols[ctr] = psi0_idx;  values[ctr++] = scal*Tdq0inv[r][2]*dI0_dpsi0;
    rows[ctr] = i_idx[r]; cols[ctr] = Eqp_idx;   values[ctr++] = scal*Tdq0inv[r][0]*dId_dEqp;
    rows[ctr] = i_idx[r]; cols[ctr] = psi1d_idx; values[ctr++] = scal*Tdq0inv[r][0]*dId_dpsi1d;
    rows[ctr] = i_idx[r]; cols[ctr] = psi2q_idx; values[ctr++] = scal*Tdq0inv[r][1]*dIq_dpsi2q;
    rows[ctr] = i_idx[r]; cols[ctr] = delta_idx;
    values[ctr++] = scal*(dTdq0inv_ddelta[r][0]*Id + dTdq0inv_ddelta[r][1]*Iq + dTdq0inv_ddelta[r][2]*I0);
    rows[ctr] = i_idx[r]; cols[ctr] = i_idx[r];  values[ctr++] = -1.0;
  }
  *nvals = ctr;
}

bool Gensal::setJacobian(gridpack::RealType **values)
{
  return true;
}

double Gensal::getInitialFieldVoltage()
{
  return Efd;
}

// Explicit mode: forward-Euler step of the rotor states with the stator
// fluxes and network voltage frozen at the start of the step.
void Gensal::preStep(double time, double timestep)
{
  if(integrationtype != EXPLICIT) return;
  double x[5], f[5];
  x[0] = Eqp; x[1] = psi1d; x[2] = psi2q; x[3] = delta; x[4] = dw;

  double theta = delta - PI/2.0;
  double tempd1 = (Xdpp - Xl)/(Xdp - Xl);
  double tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
  double param1 = (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));

  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;
  abc2dq0(vabc,time,theta,vdq0);
  double Id = (psid - tempd1*Eqp - tempd2*psi1d)/-Xdpp;
  double Iq = (psiq - psi2q)/-Xdpp;
  idq0[0] = Id; idq0[1] = Iq; idq0[2] = psi0/-Xl;

  if(hasExciter()) Efd = getExciter()->getFieldVoltage();
  if(hasGovernor()) TM = getGovernor()->getMechanicalPower();

  double dpsi1ddt = -psi1d + Eqp - (Xdp - Xl)*Id;
  LadIfd = Eqp*(1.0 + Sat(Eqp)) + (Xd - Xdp)*(Id + param1*dpsi1ddt);

  f[0] = (Efd - LadIfd)/Tdop;
  f[1] = dpsi1ddt/Tdopp;
  f[2] = (-psi2q - (Xq - Xdpp)*Iq)/Tqopp;
  f[3] = dw*OMEGA_S;
  f[4] = 1.0/(2.0*H)*((TM - D*dw)/(1 + dw) - (psid*Iq - psiq*Id));
  for(int i = 0; i < 5; i++) x[i] += timestep*f[i];
  Eqp = x[0]; psi1d = x[1]; psi2q = x[2]; delta = x[3]; dw = x[4];
}

void Gensal::postStep(double time)
{
}
