#include <regca1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Regca1::Regca1(void)
{
  nxgen   = 4; // Number of variables for this model
}

void Regca1::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgen = 4;
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
  if(pg > 0.0) {
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
  Pll_block.setparams(0.01,1.0);

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

  iabc[0] = ia*mbase/sbase;
  iabc[1] = ib*mbase/sbase;
  iabc[2] = ic*mbase/sbase;

  vabc[0] = p_va;
  vabc[1] = p_vb;
  vabc[2] = p_vc;

  Ip = Pg/p_Vm0;
  Iq = -Qg/p_Vm0;

  // Initialize blocks
  // PLL
  double Vq = Pll_block.init_given_y(OMEGA_S);

  dw = 0.0;
  theta = p_Va0;
  angle_block.init_given_y(theta);

  // Assume no limits are hit
  Ipcmd = Ip_blk.init_given_y(Ip);
  Iqcmd = Iq_blk.init_given_y(Iq);
  Vt_filter = Vt_filter_blk.init_given_u(p_Vm0);

  x[0] = p_Vm0;    // Voltage magnitude
  x[1] = pg/sbase; // Real power
  x[2] = qg/sbase; // Reactive power
  x[3] = 1.0;     // Frequency
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
    Vm  = x[0];
    Pgen  = x[1];
    Qgen  = x[2];
    Freq  = x[3];
  } else if(p_mode == XDOTVECTOBUS) {
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

  abc2dq0(vabc,time,theta,vdq0);

  Vd = vdq0[0];
  Vq = vdq0[1];

  Vt = sqrt(Vd*Vd + Vq*Vq);

  Vm_save = Vt;
  
  double omega;
  omega    = Pll_block.getoutput(Vq, timestep, true);
  dw = omega - OMEGA_S;
  theta = angle_block.getoutput(dw, timestep, true);

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

  idq0[0] = Ipout*mbase/sbase;
  idq0[1] = Iqout*mbase/sbase;
  idq0[2] = 0.0;

  getPower(time,&Pgen_save,&Qgen_save);
  Freq_save = getFreq();
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
    f[0] = Vm - Vm_save;
    f[1] = Pgen - Pgen_save;
    f[2] = Qgen - Qgen_save;
    f[3] = Freq - Freq_save;
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

  pufreq = (dw + OMEGA_S)/(2*PI*FREQ);

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

  abc2dq0(vabc,time,theta,vdq0);
  abc2dq0(iabc,time,theta,idq0);

  Vd = vdq0[0];
  Vq = vdq0[1];

  Id = idq0[0];
  Iq = idq0[1];

  V = gridpack::ComplexType(Vd,Vq);
  I = gridpack::ComplexType(Id,Iq);

  S = V*conj(I);
  Pgen = real(S);
  Qgen = imag(S);

  *Pg = Pgen*mbase/sbase;
  *Qg = Qgen*mbase/sbase;

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
  *Pg = pg/sbase;
  *Qg = qg/sbase;
}

/**
 * Return the generator current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void Regca1::getCurrent(double *ia, double *ib, double *ic)
{
  dq02abc(idq0,p_time,theta,iabc);
  
  *ia = iabc[0];
  *ib = iabc[1];
  *ic = iabc[2];
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Regca1::getCurrentGlobalLocation(int *i_gloc)
{
  if(integrationtype == IMPLICIT) {
    *i_gloc = p_gloc;
  } else {
    *i_gloc = -1; // No variables for this model
  }
}


/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values

 Number of non-zero values = 4
 */
int Regca1::matrixNumValues()
{
  int numVals;

  numVals = 4;

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
  values[ctr] = 1.0;

  rows[ctr+1] = p_gloc+1;
  cols[ctr+1] = p_gloc+1;
  values[ctr+1] = 1.0;

  rows[ctr+2] = p_gloc+2;
  cols[ctr+2] = p_gloc+2;
  values[ctr+2] = 1.0;

  rows[ctr+3] = p_gloc+3;
  cols[ctr+3] = p_gloc+3;
  values[ctr+3] = 1.0;

  ctr += 4;

  *nvals = ctr;
}
