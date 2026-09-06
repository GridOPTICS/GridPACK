/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#include <regfma1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Regfma1::Regfma1(void)
{
  nxgen = 6;
  omega = OMEGA_S;
  delta = 0.0;
  zero_Tpf = zero_Tqf = zero_Tvf = false;
  Edroop_max = Edroop_min = 0.0;
  Pmax_PI_blk_out = Pmin_PI_blk_out = Qmax_PI_blk_out = Qmin_PI_blk_out = 0.0;
  Pinv = Qinv = Vmeas = Edroop = 0.0;
  for(int i = 0; i < 3; i++) iabc[i] = diabc[i] = iout[i] = 0.0;
}

Regfma1::~Regfma1(void)
{
}

void Regfma1::getnvar(int *nvar)
{
  *nvar = nxgen;
}

void Regfma1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGenModel::load(data,idx);

  if (!data->getValue(GENERATOR_ZSOURCE, &Zsource, idx)) Zsource = gridpack::ComplexType(0.01,0.02);
  Ra = real(Zsource);
  Xl = imag(Zsource);
  L  = Xl/OMEGA_S;

  if (!data->getValue(GENERATOR_REGFMA1_MQ, &mq, idx)) mq = 0.05;
  if (!data->getValue(GENERATOR_REGFMA1_KPV, &kpv, idx)) kpv = 0.0;
  if (!data->getValue(GENERATOR_REGFMA1_KIV, &kiv, idx)) kiv = 5.86;
  if (!data->getValue(GENERATOR_REGFMA1_MP, &mp, idx)) mp = 3.77;
  if (!data->getValue(GENERATOR_REGFMA1_KPPMAX, &kppmax, idx)) kppmax = 0.05;
  if (!data->getValue(GENERATOR_REGFMA1_KIPMAX, &kipmax, idx)) kipmax = 0.2;
  if (!data->getValue(GENERATOR_REGFMA1_PMAX, &Pmax, idx)) Pmax = 1.0;
  if (!data->getValue(GENERATOR_REGFMA1_PMIN, &Pmin, idx)) Pmin = 0.0;
  if (!data->getValue(GENERATOR_REGFMA1_EMAX, &Emax, idx)) Emax = 2.0;
  if (!data->getValue(GENERATOR_REGFMA1_EMIN, &Emin, idx)) Emin = -2.0;
  if (!data->getValue(GENERATOR_REGFMA1_TPF, &Tpf, idx)) Tpf = 0.01666;
  if (!data->getValue(GENERATOR_REGFMA1_IMAX, &Imax, idx)) Imax = 2.5;
  if (!data->getValue(GENERATOR_REGFMA1_QMAX, &Qmax, idx)) Qmax = 1.0;
  if (!data->getValue(GENERATOR_REGFMA1_QMIN, &Qmin, idx)) Qmin = -1.0;
  if (!data->getValue(GENERATOR_REGFMA1_KPQMAX, &kpqmax, idx)) kpqmax = 0.1;
  if (!data->getValue(GENERATOR_REGFMA1_KIQMAX, &kiqmax, idx)) kiqmax = 10.0;
  if (!data->getValue(GENERATOR_REGFMA1_TQF, &Tqf, idx)) Tqf = 0.01666;
  if (!data->getValue(GENERATOR_REGFMA1_TVF, &Tvf, idx)) Tvf = 0.01666;
  if (!data->getValue(GENERATOR_REGFMA1_VFLAG, &Vflag, idx)) Vflag = 0;
  if (!data->getValue(GENERATOR_REGFMA1_QVFLAG, &QVflag, idx)) QVflag = 1;

  zero_Tpf = (fabs(Tpf) <= 1e-6);
  zero_Tqf = (fabs(Tqf) <= 1e-6);
  zero_Tvf = (fabs(Tvf) <= 1e-6);
  Edroop_max = Emax;
  Edroop_min = Emin;

  if(!zero_Tpf) P_filter_blk.setparams(1.0,Tpf);
  if(!zero_Tqf) Q_filter_blk.setparams(1.0,Tqf);
  if(!zero_Tvf) V_filter_blk.setparams(1.0,Tvf);
  if(Vflag == 0) Edroop_limiter_blk.setparams(1.0,Emin,Emax);
  else Edroop_PI_blk.setparams(kpv,kiv,Emin,Emax,Emin,Emax);
  Pmax_PI_blk.setparams(kppmax,kipmax,-1000.0,0.0,-1000.0,0.0);
  Pmin_PI_blk.setparams(kppmax,kipmax,0.0,1000.0,0.0,1000.0);
  Qmax_PI_blk.setparams(kpqmax,kiqmax,-1000.0,0.0,-1000.0,0.0);
  Qmin_PI_blk.setparams(kpqmax,kiqmax,0.0,1000.0,0.0,1000.0);
  Delta_blk.setparams(1.0);   // delta' = domega [rad/s]
  if(zero_Tpf || zero_Tqf || zero_Tvf)
    printf("REGFMA1 bus %d: Tpf/Tqf/Tvf = 0 samples instantaneous EMT power; results become step-size dependent\n",busnum);
}

void Regfma1::init(gridpack::RealType* xin)
{
  gridpack::RealType *x = xin + offsetb;
  double Pg = pg/mbase, Qg = qg/mbase;   // machine base
  double Vr = p_Vm0*cos(p_Va0), Vi = p_Vm0*sin(p_Va0);

  gridpack::ComplexType V(Vr,Vi), S(Pg,Qg);
  gridpack::ComplexType I = conj(S/V);
  gridpack::ComplexType E = V + Zsource*I;
  double Im = abs(I), Ia = arg(I);

  Edroop = abs(E);
  delta  = arg(E);
  omega  = OMEGA_S;
  Vmeas  = p_Vm0;
  Pinv   = Pg;
  Qinv   = Qg;

  double Vref;
  if(Vflag == 0) {
    Vref = Edroop;
  } else {
    Vref = Edroop_PI_blk.init_given_y(Edroop) + p_Vm0;
  }
  Delta_blk.init_given_y(delta);
  Pmax_PI_blk.init_given_y(0.0);
  Pmin_PI_blk.init_given_y(0.0);
  Qmax_PI_blk.init_given_y(0.0);
  Qmin_PI_blk.init_given_y(0.0);
  if(!zero_Tpf) P_filter_blk.init_given_y(Pg);
  if(!zero_Tqf) Q_filter_blk.init_given_y(Qg);
  if(!zero_Tvf) V_filter_blk.init_given_y(p_Vm0);

  // Both QVflag settings reduce to Vref = Vset + (Qset - Qinv)*mq
  if(QVflag == 0) { Vset = Vref; Qset = Qg; }
  else { Vset = Vref + Qg*mq; Qset = 0.0; }
  Pset = Pg;

  iabc[0] = Im*sin(OMEGA_S*p_time + Ia);
  iabc[1] = Im*sin(OMEGA_S*p_time + Ia - TWOPI_OVER_THREE);
  iabc[2] = Im*sin(OMEGA_S*p_time + Ia + TWOPI_OVER_THREE);
  for(int i = 0; i < 3; i++) {
    iout[i] = iabc[i]*mbase/sbase;
    x[i] = iabc[i];
    x[3+i] = iout[i];
  }
}

bool Regfma1::serialWrite(char *string, const int bufsize, const char *signal)
{
  if(!strcmp(signal,"header")) {
    sprintf(string,", %d_%s_V,%d_%s_Pg,%d_%s_Qg,%d_%s_delta, %d_%s_dw",busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str());
    return true;
  } else if(!strcmp(signal,"monitor")) {
    sprintf(string,", %6.5f,%6.5f,%6.5f,%6.5f, %6.5f",Vmeas,Pinv*mbase/sbase,Qinv*mbase/sbase,delta,getSpeedDeviation());
    return true;
  }
  return false;
}

void Regfma1::write(const char* signal, char* string)
{
}

void Regfma1::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values + offsetb;
  if(p_mode == XVECTOBUS) {
    for(int i = 0; i < 3; i++) { iabc[i] = x[i]; iout[i] = x[3+i]; }
  } else if(p_mode == XDOTVECTOBUS) {
    for(int i = 0; i < 3; i++) diabc[i] = x[i];
  }
}

void Regfma1::preStep(double time, double timestep)
{
  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;
  abc2dq0(vabc,time,delta,vdq0);
  abc2dq0(iabc,time,delta,idq0);

  double Vd = vdq0[0], Vq = vdq0[1], Id = idq0[0], Iq = idq0[1];
  double Vt = sqrt(Vd*Vd + Vq*Vq);
  double Pgen = Vd*Id + Vq*Iq;
  double Qgen = Vq*Id - Vd*Iq;

  Pinv  = zero_Tpf ? Pgen : P_filter_blk.getoutput(Pgen,timestep,true);
  Qinv  = zero_Tqf ? Qgen : Q_filter_blk.getoutput(Qgen,timestep,true);
  Vmeas = zero_Tvf ? Vt   : V_filter_blk.getoutput(Vt,timestep,true);

  // Current limit: cap E so that |I| cannot exceed Imax (same as phasor model)
  double Im = sqrt(Id*Id + Iq*Iq);
  if(Im >= Imax) {
    gridpack::ComplexType V(Vd,Vq), Ilim(Imax*Id/Im, Imax*Iq/Im);
    Edroop_max = abs(V + Zsource*Ilim);
    Edroop_min = 0.0;
  } else {
    Edroop_max = Emax;
    Edroop_min = Emin;
  }

  Qmax_PI_blk_out = Qmax_PI_blk.getoutput(Qmax - Qinv,timestep,true);
  Qmin_PI_blk_out = Qmin_PI_blk.getoutput(Qmin - Qinv,timestep,true);
  double Vref = Vset + (Qset - Qinv)*mq + Qmax_PI_blk_out + Qmin_PI_blk_out;
  if(Vflag == 0) {
    Edroop = Edroop_limiter_blk.getoutput(Vref,Edroop_min,Edroop_max);
  } else {
    Edroop = Edroop_PI_blk.getoutput(Vref - Vmeas,timestep,-1000.0,1000.0,Edroop_min,Edroop_max,true);
  }

  Pmax_PI_blk_out = Pmax_PI_blk.getoutput(Pmax - Pinv,timestep,true);
  Pmin_PI_blk_out = Pmin_PI_blk.getoutput(Pmin - Pinv,timestep,true);
  double domega = OMEGA_S*(mp*(Pset - Pinv) + Pmax_PI_blk_out + Pmin_PI_blk_out);
  omega = OMEGA_S + domega;
  delta = Delta_blk.getoutput(domega,timestep,true);
}

void Regfma1::postStep(double time)
{
}

void Regfma1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values + offsetb;
  if(p_mode != RESIDUAL_EVAL) return;

  eabc[0] = Edroop*sin(OMEGA_S*p_time + delta);
  eabc[1] = Edroop*sin(OMEGA_S*p_time + delta - TWOPI_OVER_THREE);
  eabc[2] = Edroop*sin(OMEGA_S*p_time + delta + TWOPI_OVER_THREE);
  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;
  for(int i = 0; i < 3; i++) {
    f[i]   = eabc[i] - vabc[i] - Ra*iabc[i] - L*diabc[i];
    f[3+i] = iabc[i]*mbase/sbase - iout[i];
  }
}

double Regfma1::getFreq()
{
  return omega/OMEGA_S;
}

void Regfma1::getPower(double time, double *Pgen, double *Qgen)
{
  *Pgen = Pinv*mbase/sbase;
  *Qgen = Qinv*mbase/sbase;
}

void Regfma1::getInitialPower(double *Pg, double *Qg)
{
  *Pg = pg/sbase;
  *Qg = qg/sbase;
}

void Regfma1::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = iout[0]; *ib = iout[1]; *ic = iout[2];
}

void Regfma1::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc + 3;
}

int Regfma1::matrixNumValues()
{
  return 12;
}

void Regfma1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  for(int i = 0; i < 3; i++) {
    rows[ctr] = p_gloc + i; cols[ctr] = p_gloc + i;        values[ctr++] = -Ra - L*shift;
    rows[ctr] = p_gloc + i; cols[ctr] = p_glocvoltage + i; values[ctr++] = -1.0;
  }
  for(int i = 0; i < 3; i++) {
    rows[ctr] = p_gloc + 3 + i; cols[ctr] = p_gloc + i;     values[ctr++] = mbase/sbase;
    rows[ctr] = p_gloc + 3 + i; cols[ctr] = p_gloc + 3 + i; values[ctr++] = -1.0;
  }
  *nvals = ctr;
}
