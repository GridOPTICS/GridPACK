#include <gdform.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Gdform::Gdform(void)
{
  nxgen   = 6; // Number of variables for this model
}

void Gdform::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgen = 6;
  *nvar = nxgen;
}

Gdform::~Gdform(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseGeneratorModel
 */
void Gdform::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGenModel::load(data,idx); // load parameters in base generator model

  // load parameters for the model type
  data->getValue(GENERATOR_ZSOURCE,&Zsource,idx);
  L = imag(Zsource)/OMEGA_S;

    if (!data->getValue(GENERATOR_MQ, &mq, idx)) mq=0.05; // 
  if (!data->getValue(GENERATOR_KPV, &kpv, idx)) kpv=0.0; // 
  if (!data->getValue(GENERATOR_KIV, &kiv, idx)) kiv=5.86; // 
  if (!data->getValue(GENERATOR_MP, &mp, idx)) mp=3.77; // 
  if (!data->getValue(GENERATOR_KPPMAX, &kppmax, idx)) kppmax=0.05; // 
  if (!data->getValue(GENERATOR_KIPMAX, &kipmax, idx)) kipmax=0.2; // 
  if (!data->getValue(GENERATOR_PMAX, &Pmax, idx)) Pmax=1.0; // 
  if (!data->getValue(GENERATOR_PMIN, &Pmin, idx)) Pmin=0.0; // 

  if (!data->getValue(GENERATOR_EMAX, &Emax, idx)) Emax = 2.0; // 
  if (!data->getValue(GENERATOR_EMIN, &Emin, idx)) Emin = -2.0; //

  Edroop_min = Emin;
  Edroop_max = Emax;

  if (!data->getValue(GENERATOR_TPF, &Tpf, idx)) Tpf = 0.01666; //
  if(fabs(Tpf) < 1e-6) zero_Tpf = true;

  if (!data->getValue(GENERATOR_IMAX, &Imax, idx)) Imax=2.5;
  if (!data->getValue(GENERATOR_QMAX, &Qmax, idx)) Qmax=1.0;
  if (!data->getValue(GENERATOR_QMIN, &Qmin, idx)) Qmin=-1.0;

  if (!data->getValue(GENERATOR_KPQMAX, &kpqmax, idx)) kpqmax=0.1; // 
  if (!data->getValue(GENERATOR_KIQMAX, &kiqmax, idx)) kiqmax=10; // 

  if (!data->getValue(GENERATOR_TQF, &Tqf, idx)) Tqf=0.01666; //
  if(fabs(Tqf) < 1e-6) zero_Tqf = true;

  if (!data->getValue(GENERATOR_TVF, &Tvf, idx)) Tvf=0.01666; // 
  if(fabs(Tvf) < 1e-6) zero_Tvf = true;

  if(!data->getValue(GENERATOR_VFLAG,&Vflag, idx)) Vflag = 0;

  // Initialize blocks
  if(!zero_Tpf) {
    P_filter_blk.setparams(1.0,Tpf);
  }

  if(!zero_Tqf) {
    Q_filter_blk.setparams(1.0,Tqf);
  }

  if(!zero_Tvf) {
    V_filter_blk.setparams(1.0,Tvf);
  }

  if(Vflag == 0) {
    Edroop_limiter_blk.setparams(1.0,Emin,Emax);
  } else {
    Edroop_PI_blk.setparams(kpv,kiv);
  }

  Pmax_PI_blk.setparams(kppmax,kipmax,-1000.0,0.0,-1000.0,0.0);
  Pmin_PI_blk.setparams(kppmax,kipmax,0.0,1000.0,0.0,1000.0);

  Qmax_PI_blk.setparams(kpqmax,kiqmax,-1000.0,0.0,-1000.0,0.0);
  Qmin_PI_blk.setparams(kpqmax,kiqmax,0.0,1000.0,0.0,1000.0);

  Delta_blk.setparams(1.0);
}

/**
 * Initialize generator model before calculation
 * @param [output] xin - array where initialized generator variables should be set
 */
void Gdform::init(gridpack::RealType* xin)
{
  double Pg, Qg;  // Generator real and reactive power
  gridpack::RealType *x = xin+offsetb; // generator array starts from this location
  double Vr, Vi;

  Pg = pg/mbase;
  Qg = qg/mbase;

  Vr = p_Vm0*cos(p_Va0);
  Vi = p_Vm0*sin(p_Va0);

  gridpack::ComplexType V = gridpack::ComplexType(Vr,Vi);
  gridpack::ComplexType S = gridpack::ComplexType(Pg,Qg);
  gridpack::ComplexType I;

  // Machine output current
  I = conj(S/V);
  double Im = abs(I);
  double Ia = arg(I);

  E = V + Zsource*I;

  Edroop = abs(E);
  delta  = arg(E);

  // Initialize blocks
  double Vref;
  if(Vflag == 0) {
    Vref = Edroop;
  } else {
    Vref = Edroop_PI_blk.init_given_y(Edroop);
    Vref += p_Vm0; 
  }
  
  Delta_blk.init_given_y(delta);

  Pmax_PI_blk.init_given_y(0.0);
  Pmin_PI_blk.init_given_y(0.0);

  Qmax_PI_blk.init_given_y(0.0);
  Qmin_PI_blk.init_given_y(0.0);

  if(!zero_Tpf) {
    P_filter_blk.init_given_y(Pg);
  }

  if(!zero_Tqf) {
    Q_filter_blk.init_given_y(Qg);
  }

  if(!zero_Tvf) {
    V_filter_blk.init_given_y(p_Vm0);
  }
  
  Vset = Vref + Qg*mq;
  Pset = Pg;

  double ia,ib,ic;
  ia = Im*sin(OMEGA_S*p_time + Ia);
  ib = Im*sin(OMEGA_S*p_time + Ia - TWOPI_OVER_THREE);
  ic = Im*sin(OMEGA_S*p_time + Ia + TWOPI_OVER_THREE);

  iabc[0] = ia;
  iabc[1] = ib;
  iabc[2] = ic;

  vabc[0] = p_va;
  vabc[1] = p_vb;
  vabc[2] = p_vc;

  x[0] = ia;
  x[1] = ib;
  x[2] = ic;
  x[3] = iabc[0]*mbase/sbase; // current on system MVAbase
  x[4] = iabc[1]*mbase/sbase; // current on system MVAbase
  x[5] = iabc[2]*mbase/sbase; // current on system MVAbase
}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Gdform::serialWrite(char *string, const int bufsize,const char *signal)
{
  if(!strcmp(signal,"header")) {
    /* Print output header */
    sprintf(string,", %d_%s_V,%d_%s_Pg,%d_%s_delta, %d_%s_dw",busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str());
    return true;
  } else if(!strcmp(signal,"monitor")) {
    /* Print output */
    sprintf(string,", %6.5f,%6.5f,%6.5f, %6.5f",Vmeas,Pinv*mbase/sbase,delta,(omega-OMEGA_S)/OMEGA_S);
    return true;
  }
  return false;
}


/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Gdform::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void Gdform::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values+offsetb; // generator array starts from this location

  if(p_mode == XVECTOBUS) {
    iabc[0] = x[0];
    iabc[1] = x[1];
    iabc[2] = x[2];
    iout[0] = x[3];
    iout[1] = x[4];
    iout[2] = x[5];
  } else if(p_mode == XDOTVECTOBUS) {
    diabc[0] = x[0];
    diabc[1] = x[1];
    diabc[2] = x[2];
  }
  
}

/**
   Prestep function
*/
void Gdform::preStep(double time ,double timestep)
{
  double Vd,Vq, Vt;

  vabc[0] = p_va;
  vabc[1] = p_vb;
  vabc[2] = p_vc;

  abc2dq0(vabc,time,delta,vdq0);
  abc2dq0(iabc,time,delta,idq0);

  Vd = vdq0[0];
  Vq = vdq0[1];

  Vt = sqrt(Vd*Vd + Vq*Vq);

  double Pgen, Qgen;
  Pgen = vdq0[0]*idq0[0] + vdq0[1]*idq0[1];
  Qgen = vdq0[1]*idq0[0] - vdq0[0]*idq0[1];

  if(!zero_Tpf) {
    Pinv = P_filter_blk.getoutput(Pgen,timestep,true);
  } else {
    Pinv = Pgen;
  }

  if(!zero_Tqf) {
    Qinv = Q_filter_blk.getoutput(Qgen,timestep,true);
  } else {
    Qinv = Qgen;
  }

  if(!zero_Tvf) {
    Vmeas = V_filter_blk.getoutput(Vt, timestep,true);
  } else {
    Vmeas = Vt;
  }

  Qmax_PI_blk_out = Qmax_PI_blk.getoutput(Qmax-Qinv,timestep,true);
  Qmin_PI_blk_out = Qmin_PI_blk.getoutput(Qmin-Qinv,timestep,true);

  double Vref;
  Vref = Vset - Qinv*mq + Qmax_PI_blk_out + Qmin_PI_blk_out;
  if(Vflag == 0) {
    Edroop = Edroop_limiter_blk.getoutput(Vref,Edroop_min,Edroop_max);
  } else {
    Edroop = Edroop_PI_blk.getoutput(Vref-Vmeas,timestep,-1000.0,1000.0,Edroop_min,Edroop_max,true);
  }

  Pmax_PI_blk_out = Pmax_PI_blk.getoutput(Pmax-Pinv,timestep,true);
  Pmin_PI_blk_out = Pmin_PI_blk.getoutput(Pmin-Pinv,timestep,true);

  
  double domega;

  domega = OMEGA_S*(mp*(Pset - Pinv) + Pmax_PI_blk_out + Pmin_PI_blk_out);

  omega = OMEGA_S + domega;

  delta = Delta_blk.getoutput(domega,timestep,true);

}

/**
   Poststep function
*/
void Gdform::postStep(double time)
{
}



/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Gdform::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // generator array starts from this location

  if(p_mode == RESIDUAL_EVAL) {
    eabc[0] = Edroop*sin(OMEGA_S*p_time + delta);
    eabc[1] = Edroop*sin(OMEGA_S*p_time + delta - 2.0*PI/3.0);
    eabc[2] = Edroop*sin(OMEGA_S*p_time + delta + 2.0*PI/3.0);

    vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;
    
    f[0] = eabc[0] - vabc[0] - L*diabc[0];
    f[1] = eabc[1] - vabc[1] - L*diabc[1];
    f[2] = eabc[2] - vabc[2] - L*diabc[2];

    f[3] = iabc[0]*mbase/sbase - iout[0];
    f[4] = iabc[1]*mbase/sbase - iout[1];
    f[5] = iabc[2]*mbase/sbase - iout[2];
  }
}

/**
 * Return the generator frequency (pu)
 * @param [output] freq - machine frequency
 *
 * Note: Frequency is per unit. Steady-state frequency is 1.0
 */
double Gdform::getFreq()
{
  double pufreq;

  return pufreq;
}


/**
 * Return the generator real and reactive power
 * @param [input] time - the current time
 * @param [output] Pg - generator real power
 * @param [output] Qg - generator reactive power
 *
 * Note: Power is on system MVA base
 */
void Gdform::getPower(double time,double *Pgen, double *Qgen)
{

  *Pgen = Pinv*mbase/sbase;
  *Qgen = Qinv*mbase/sbase;

}

/**
 * Return the generator initial real and reactive power
 * @param [output] Pg(t0) - generator real power
 * @param [output] Qg(t0) - generator reactive power
 *
 * Note: Power is pu on system MVA base
 */
void Gdform::getInitialPower(double *Pg, double *Qg)
{
  *Pg = pg/sbase;
  *Qg = qg/sbase;
}

/**
 * Return the generator current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void Gdform::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = iout[0];
  *ib = iout[1];
  *ic = iout[2];
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Gdform::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc+3;
}


/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values

 Number of non-zero values = 4
 */
int Gdform::matrixNumValues()
{
  int numVals;

  numVals = 12;

  return numVals;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Gdform::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  rows[ctr] = p_gloc;
  cols[ctr] = p_gloc;
  values[ctr] = -L*shift;

  rows[ctr+1] = p_gloc;
  cols[ctr+1] = p_glocvoltage;
  values[ctr+1] = -1.0;

  ctr += 2;

  rows[ctr] = p_gloc+1;
  cols[ctr] = p_gloc+1;
  values[ctr] = -L*shift;

  rows[ctr+1] = p_gloc+1;
  cols[ctr+1] = p_glocvoltage+1;
  values[ctr+1] = -1.0;

  ctr += 2;

  rows[ctr] = p_gloc+2;
  cols[ctr] = p_gloc+2;
  values[ctr] = -L*shift;

  rows[ctr+1] = p_gloc+2;
  cols[ctr+1] = p_glocvoltage+2;
  values[ctr+1] = -1.0;

  ctr += 2;

  rows[ctr] = p_gloc + 3;
  cols[ctr] = p_gloc;
  values[ctr] = mbase/sbase;

  rows[ctr+1] = p_gloc + 3;
  cols[ctr+1] = p_gloc + 3;
  values[ctr+1] = -1.0;

  ctr += 2;

  rows[ctr] = p_gloc + 4;
  cols[ctr] = p_gloc + 1;
  values[ctr] = mbase/sbase;

  rows[ctr+1] = p_gloc + 4;
  cols[ctr+1] = p_gloc + 4;
  values[ctr+1] = -1.0;

  ctr += 2;

  rows[ctr] = p_gloc + 5;
  cols[ctr] = p_gloc + 2;
  values[ctr] = mbase/sbase;

  rows[ctr+1] = p_gloc + 5;
  cols[ctr+1] = p_gloc + 5;
  values[ctr+1] = -1.0;

  ctr += 2;

  *nvals = ctr;
}
