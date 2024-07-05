#include <regca1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Regca1::Regca1(void)
{
  nxgen   = 6; // Number of variables for this model
}

void Regca1::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgen = 6;
  *nvar = nxgen;
}

Regca1::~Regca1(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseGeneratorModel
 */
void Regca1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGenModel::load(data,idx); // load parameters in base generator model

  // load parameters for the model type
  data->getValue(GENERATOR_ZSOURCE,&Zsource,idx);
  L = imag(Zsource)/OMEGA_S;

  if (!data->getValue(GENERATOR_REGC_LVPLSW , &lvplsw, idx)) lvplsw = 0; 
  if (!data->getValue(GENERATOR_REGC_TG  ,    &tg, idx))     tg = 0.02; 
  if (!data->getValue(GENERATOR_REGC_RRPWR  , &rrpwr, idx))  rrpwr = 10.0; 
  if (!data->getValue(GENERATOR_REGC_BRKPT  , &brkpt, idx))  brkpt = 0.9; 
  if (!data->getValue(GENERATOR_REGC_ZEROX  , &zerox, idx))  zerox = 0.4; 
  if (!data->getValue(GENERATOR_REGC_LVPL1  , &lvpl1, idx))  lvpl1 = 1.22; 
  if (!data->getValue(GENERATOR_REGC_VOLIM  , &volim, idx))  volim = 1.2; 
  if (!data->getValue(GENERATOR_REGC_LVPNT1 , &lvpnt1, idx)) lvpnt1 = 0.8; 
  if (!data->getValue(GENERATOR_REGC_LVPNT0 , &lvpnt0, idx)) lvpnt0 = 0.4; 
  if (!data->getValue(GENERATOR_REGC_LOLIM  , &lolim, idx))  lolim = -1.3; 
  if (!data->getValue(GENERATOR_REGC_TFLTR  , &tfltr, idx))  tfltr = 0.02; 
  if (!data->getValue(GENERATOR_REGC_KHV  ,   &khv, idx))    khv = 0.0; 
  if (!data->getValue(GENERATOR_REGC_IQRMAX , &iqrmax, idx)) iqrmax = 999.0; 
  if (!data->getValue(GENERATOR_REGC_IQRMIN , &iqrmin, idx)) iqrmin = -999.0; 
  if (!data->getValue(GENERATOR_REGC_ACCEL  , &accel, idx))  accel = 0.7;

  // Set up blocks

  // transfer function blocks
  Ip_blk.setparams(1.0,tg);
  
  Iq_blk.setparams(-1.0,tg);
  if(qg > 0.0) {
    // Upper limit active with Qg > 0
    Iq_blk.setdxlimits(-1000.0,iqrmax);
  } else {
    // Lower limit active when Qg < 0
    Iq_blk.setdxlimits(iqrmin,1000.0);
  }
  Vt_filter_blk.setparams(1.0,tfltr);

  double u[2],y[2];

  // Lvpl block
  u[0] = zerox; u[1] = brkpt;
  y[0] = 0.0;   y[1] = lvpl1;

  Lvpl_blk.setparams(u,y,y[0],1000.0); // Infinite upper limit

  // Lvpnt block
  u[0] = lvpnt0; u[1] = lvpnt1;
  y[0] = 0.0;    y[1] = 1.0;

  Lvpnt_blk.setparams(u,y,y[0],y[1]);

  // Iq low limiter
  Iqlowlim_blk.setparams(1.0,lolim,1000.0);

  // PLL block
  omega_Pll_block.setparams(0.1,1.0);

  // Integrator block
  angle_block.setparams(1.0);

}

/**
 * Initialize generator model before calculation
 * @param [output] xin - array where initialized generator variables should be set
 */
void Regca1::init(gridpack::RealType* xin)
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

  Ip = Pg/p_Vm0;
  Iq = -Qg/p_Vm0;

  Ipref = Ip;
  Iqref = Iq;

  // Initialize blocks

  domega = 0.0;
  // PLL
  double Vq = omega_Pll_block.init_given_y(domega);

  omega = OMEGA_S*(1 + domega);
  
  delta = delta0 = p_Va0;
  angle_block.init_given_y(delta);

  // Assume no limits are hit
  Ipcmd = Ip_blk.init_given_y(Ip);
  Iqcmd = Iq_blk.init_given_y(Iq);
  Vt_filter = Vt_filter_blk.init_given_u(p_Vm0);

  x[0] = ia;
  x[1] = ib;
  x[2] = ic;
  x[3] = ia*mbase/sbase; // current on system base
  x[4] = ib*mbase/sbase; // current on system base
  x[5] = ic*mbase/sbase; // current on system base
}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Regca1::serialWrite(char *string, const int bufsize,const char *signal)
{
  if(!strcmp(signal,"header")) {
    /* Print output header */
    sprintf(string,", %d_%s_V,%d_%s_Pg,%d_%s_delta, %d_%s_dw",busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str());
    return true;
  } else if(!strcmp(signal,"monitor")) {
    /* Print output */
    getPower(p_time,&Pgen,&Qgen);
    sprintf(string,", %6.5f,%6.5f,%6.5f, %6.5f",Vt_filter,Pgen,delta,domega);
    return true;
  }
  return false;
}


/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Regca1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void Regca1::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values+offsetb; // generator array starts from this location

  if(p_mode == XVECTOBUS) {
    iabc[0]  = x[0];
    iabc[1]  = x[1];
    iabc[2]  = x[2];
    iout[0]  = x[3];
    iout[1]  = x[4];
    iout[2]  = x[5];
  } else if(p_mode == XDOTVECTOBUS) {
    diabc[0] = x[0];
    diabc[1] = x[1];
    diabc[2] = x[2];
  }
  
}

/**
   Prestep function
*/
void Regca1::preStep(double time ,double timestep)
{
  double Vd,Vq, Vt;
  double Iq_olim;

  vabc[0] = p_va;
  vabc[1] = p_vb;
  vabc[2] = p_vc;

  abc2dq0(vabc,time,delta,vdq0);

  Vd = vdq0[0];
  Vq = vdq0[1];

  Vt = sqrt(Vd*Vd + Vq*Vq);

  gridpack::ComplexType Vdq = gridpack::ComplexType(Vd,Vq);
  double ang;
  ang = arg(Vdq);
  domega = omega_Pll_block.getoutput(ang, timestep, true);
  //  omega  = OMEGA_S*(1 + domega);
  delta  = angle_block.getoutput(OMEGA_S*domega, timestep, true);

  Vt_filter = Vt_filter_blk.getoutput(Vt, timestep, true);

  if(hasExciter()) {
    getExciter()->getIpcmdIqcmd(&Ipcmd,&Iqcmd);
  }
  
  if(lvplsw) {
    Lvpl_out = Lvpl_blk.getoutput(Vt_filter);

    Ip = Ip_blk.getoutput(Ipcmd,timestep,-1000.0,Lvpl_out,-1000.0,rrpwr,-1000.0,1000.0,true);
  } else {
    Ip = Ip_blk.getoutput(Ipcmd,timestep,-1000.0,1000.0,-1000.0,rrpwr,-1000.0,1000.0,true);
  }
				   
  Iq = Iq_blk.getoutput(Iqcmd,timestep,true);

  Lvpnt_out = Lvpnt_blk.getoutput(Vt);

  Ipout = Ip*Lvpnt_out;

  Iq_olim = std::max(0.0,khv*(Vt - volim));

  Iqout = Iqlowlim_blk.getoutput(Iq - Iq_olim);
}

/**
   Poststep function
*/
void Regca1::postStep(double time)
{
}



/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Regca1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // generator array starts from this location

  if(p_mode == RESIDUAL_EVAL) {
    double Ipq0[3], igen[3];
    Ipq0[0] = Ip;
    Ipq0[1] = Iq;
    Ipq0[2] = 0.0;

    dq02abc(Ipq0,p_time, delta, igen);

    f[0] = igen[0] - iabc[0];
    f[1] = igen[1] - iabc[1];
    f[2] = igen[2] - iabc[2];
    
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
double Regca1::getFreq()
{
  double pufreq;

  pufreq = 1 + domega;

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
void Regca1::getPower(double time,double *Pg, double *Qg)
{
  double Vd, Vq, Id, Iq;
  gridpack::ComplexType V,I,S;
  double idq0[3];
  
  vabc[0] = p_va;
  vabc[1] = p_vb;
  vabc[2] = p_vc;

  abc2dq0(vabc,time,delta,vdq0);
  abc2dq0(iabc,time,delta,idq0);

  Vd = vdq0[0];
  Vq = vdq0[1];

  Id = idq0[0];
  Iq = idq0[1];

  V = gridpack::ComplexType(Vd,Vq);
  I = gridpack::ComplexType(Id,Iq);

  S = V*conj(I);
  Pgen = real(S);
  Qgen = imag(S);

  *Pg = Pgen;
  *Qg = Qgen;

}

/**
 * Return the generator initial real and reactive power
 * @param [output] Pg(t0) - generator real power
 * @param [output] Qg(t0) - generator reactive power
 *
 * Note: Power is pu on system MVA base
 */
void Regca1::getInitialPower(double *Pg, double *Qg)
{
  *Pg = pg/mbase;
  *Qg = qg/mbase;
}

/**
 * Return the generator current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void Regca1::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = iout[0];
  *ib = iout[1];
  *ic = iout[2];
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Regca1::getCurrentGlobalLocation(int *i_gloc)
{
    *i_gloc = p_gloc + 3;
}


/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values

 Number of non-zero values = 12
 */
int Regca1::matrixNumValues()
{
  int numVals;

  numVals = 9;

  return numVals;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Regca1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  rows[ctr] = p_gloc;
  cols[ctr] = p_gloc;
  values[ctr] = -1.0;

  rows[ctr+1] = p_gloc+1;
  cols[ctr+1] = p_gloc+1;
  values[ctr+1] = -1.0;

  rows[ctr+2] = p_gloc+2;
  cols[ctr+2] = p_gloc+2;
  values[ctr+2] = -1.0;

  ctr += 3;

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
