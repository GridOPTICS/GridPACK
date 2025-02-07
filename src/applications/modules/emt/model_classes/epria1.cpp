#include <epria1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Epria1::Epria1(void)
{
  nxgen   = 6; // Number of variables for this model when integration type is implicit or implicit-explicit
}

void Epria1::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgen = 6;
  *nvar = nxgen;
}

Epria1::~Epria1(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseGeneratorModel
 */
void Epria1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGenModel::load(data,idx); // load parameters in base generator model

  Model_PrintInfo();

  // load parameters for the model type

  if(!data->getValue(EPRIA1_PARAM1, &modelparams.Vbase, idx)) modelparams.Vbase = 0.0;
  if(!data->getValue(EPRIA1_PARAM2, &modelparams.Sbase, idx)) modelparams.Sbase = 0.0;
  if(!data->getValue(EPRIA1_PARAM3, &modelparams.Vdcbase, idx)) modelparams.Vdcbase = 0.0;
  if(!data->getValue(EPRIA1_PARAM4, &modelparams.KpI, idx)) modelparams.KpI = 0.0;
  if(!data->getValue(EPRIA1_PARAM5, &modelparams.KiI, idx)) modelparams.KiI = 0.0;
  if(!data->getValue(EPRIA1_PARAM6, &modelparams.KpPLL, idx)) modelparams.KpPLL = 0.0;
  if(!data->getValue(EPRIA1_PARAM7, &modelparams.KiPLL, idx)) modelparams.KiPLL = 0.0;
  if(!data->getValue(EPRIA1_PARAM8, &modelparams.KpP, idx)) modelparams.KpP = 0.0;
  if(!data->getValue(EPRIA1_PARAM9, &modelparams.KiP, idx)) modelparams.KiP = 0.0;
  if(!data->getValue(EPRIA1_PARAM10, &modelparams.KpQ, idx)) modelparams.KpQ = 0.0;
  if(!data->getValue(EPRIA1_PARAM11, &modelparams.KiQ, idx)) modelparams.KiQ = 0.0;
  if(!data->getValue(EPRIA1_PARAM12, &modelparams.Imax, idx)) modelparams.Imax = 0.0;
  if(!data->getValue(EPRIA1_PARAM13, &modelparams.PQflag, idx)) modelparams.PQflag = 0;
  if(!data->getValue(EPRIA1_PARAM14, &modelparams.Vdip, idx)) modelparams.Vdip = 0.0;
  if(!data->getValue(EPRIA1_PARAM15, &modelparams.Vup, idx)) modelparams.Vup = 0.0;
  if(!data->getValue(EPRIA1_PARAM16, &modelparams.Rchoke, idx)) modelparams.Rchoke = 0.0;
  if(!data->getValue(EPRIA1_PARAM17, &modelparams.Lchoke, idx)) modelparams.Lchoke = 0.0;
  if(!data->getValue(EPRIA1_PARAM18, &modelparams.Cfilt, idx)) modelparams.Cfilt = 0.0;
  if(!data->getValue(EPRIA1_PARAM19, &modelparams.Rdamp, idx)) modelparams.Rdamp = 0.0;

  zero_Rdamp = false;
  if(modelparams.Rdamp < 1e-6) zero_Rdamp = true;
  

  // Overwrite Zsource with values from the model
  Zsource = gridpack::ComplexType(modelparams.Rchoke,modelparams.Lchoke);
  
  model.Parameters      = (void*)&modelparams;
  model.ExternalInputs  = (void*)&modelinputs;
  model.ExternalOutputs = (void*)&modeloutputs;
  model.DoubleStates    = (double*)modelstates;

  int ok = Model_CheckParameters(&model);

  p_Rs = modelparams.Rchoke;
  p_L  = modelparams.Lchoke/OMEGA_S;
}

/**
 * Initialize generator model before calculation
 * @param [output] xin - array where initialized generator variables should be set
 */
void Epria1::init(gridpack::RealType* xin)
{
  double IGD,IGQ; // Machine currents in cartesian coordinates
  double Pg, Qg;  // Generator real and reactive power
  double dw=0.0;  // Initial machine speed deviation
  gridpack::RealType *x = xin+offsetb; // generator array starts from this location
  double Rdamp = modelparams.Rdamp;
  double Cfilt = modelparams.Cfilt;

  Pg = pg/mbase;
  Qg = qg/mbase;

  VD = p_Vm0*cos(p_Va0);
  VQ = p_Vm0*sin(p_Va0);

  gridpack::ComplexType V = gridpack::ComplexType(VD,VQ);
  gridpack::ComplexType S = gridpack::ComplexType(Pg,Qg);
  gridpack::ComplexType I,I_RC,IL,IC;
  gridpack::ComplexType E;
  gridpack::ComplexType j = gridpack::ComplexType(0.0,1.0);

  modelstates[2] = p_Va0;

  modelinputs.Va = p_va;
  modelinputs.Vb = p_vb;
  modelinputs.Vc = p_vc;

  modelinputs.Pref = Pg;
  modelinputs.Qref = Qg;

  modelinputs.Vref = abs(V);

  // Total current
  I = conj(S/V);

  // Current through RC filter
  I_RC = (zero_Rdamp?0.0:V/Rdamp) - j*(V*Cfilt);

  IC = j*(V*Cfilt);
  IL = I - (zero_Rdamp?0.0:V/Rdamp) + IC;
  
  E = V + IL*Zsource;

  double Im, Ia;
  double ILm,ILa;
  double ICm,ICang;
  double ia,ib,ic;
  double iLa,iLb,iLc;
  double iCa,iCb,iCc;
  double iRa,iRb,iRc;
  double ea,eb,ec;
  double Em,Ea;
  double e[3];

  Em = abs(E);
  Ea = arg(E);

  ea = e[0] = Em*sin(Ea);
  eb = e[1] = Em*sin(Ea - 2*PI/3.0);
  ec = e[2] = Em*sin(Ea + 2*PI/3.0);

  Im = abs(I);
  Ia = arg(I);
  ILm = abs(IL);
  ILa = arg(IL);
  ICm   = abs(IC);
  ICang = arg(IC);

  ia = Im*sin(Ia);
  ib = Im*sin(Ia - 2*PI/3.0);
  ic = Im*sin(Ia + 2*PI/3.0);

  iLa = ILm*sin(ILa);
  iLb = ILm*sin(ILa - 2*PI/3.0);
  iLc = ILm*sin(ILa + 2*PI/3.0);

  iCa = ICm*sin(ICang);
  iCb = ICm*sin(ICang - 2*PI/3.0);
  iCc = ICm*sin(ICang + 2*PI/3.0);

  iRa = (zero_Rdamp?0.0:p_va/Rdamp);
  iRb = (zero_Rdamp?0.0:p_vb/Rdamp);
  iRc = (zero_Rdamp?0.0:p_vc/Rdamp);

  p_vabc[0] = p_va;
  p_vabc[1] = p_vb;
  p_vabc[2] = p_vc;

  double i[3];
  i[0] = ia;
  i[1] = ib;
  i[2] = ic;

  abc2dq0(p_vabc,p_time,p_Va0,p_vdq0);
  abc2dq0(i,p_time,p_Va0,p_idq0);

  double Edq0[3];

  abc2dq0(e,p_time,p_Va0,Edq0);
  Ed = Edq0[0];
  Eq = Edq0[1];
  
  // Idref and Iqref
  double iLdq0ref[3],iL[3];

  iL[0] = iLa;
  iL[1] = iLb;
  iL[2] = iLc;

  abc2dq0(iL,p_time,p_Va0,iLdq0ref);

  // Expected outputs at t = 0
  modeloutputs.Idrefout = p_idq0[0];
  modeloutputs.Iqrefout = p_idq0[1];
  modeloutputs.Freqpll = 60.0;

  x[0] = p_iabc[0] = ia*mbase/sbase;
  x[1] = p_iabc[1] = ib*mbase/sbase;
  x[2] = p_iabc[2] = ic*mbase/sbase;
  x[3] = p_iLabc[0] = iLa;
  x[4] = p_iLabc[1] = iLb;
  x[5] = p_iLabc[2] = iLc;

  modeloutputs.Ea = ea;
  modeloutputs.Eb = eb;
  modeloutputs.Ec = ec;

  modeloutputs.Pout = Pg;
  modeloutputs.Qout = Qg;

  modelinputs.Ia = ia;
  modelinputs.Ib = ib;
  modelinputs.Ic = ic;

  modelinputs.IaL1 = iLa;
  modelinputs.IbL1 = iLb;
  modelinputs.IcL1 = iLc;

  double Cshunt[3][3];
  Cshunt[0][0] = Cshunt[1][1] = Cshunt[2][2] = Cfilt/OMEGA_S;
  bus->addCshunt(Cshunt,1.0);

  int ok = Model_Initialize(&model);

}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Epria1::serialWrite(char *string, const int bufsize,const char *signal)
{
  if(!strcmp(signal,"header")) {
    /* Print output header */
    sprintf(string,", %d_%s_V,%d_%s_Pg,%d_%s_delta, %d_%s_dw",busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str());
    return true;
  } else if(!strcmp(signal,"monitor")) {
    /* Print output */
      double Vm, Pgen,dspd=0.0,vdq0[3],idq0[3];
      abc2dq0(p_vabc,p_time,phi_PLL,vdq0);
      abc2dq0(p_iabc,p_time,phi_PLL,idq0);
    if(p_online) {
      Vm = sqrt(vdq0[0]*vdq0[0] + vdq0[1]*vdq0[1]);
    } else {
      Vm = p_Vm0;
    }
    Pgen = p_online*(vdq0[0]*idq0[0] + vdq0[1]*idq0[1]);
    sprintf(string,", %6.5f,%6.5f,%6.5f, %6.5f",Vm,Pgen,phi_PLL,modeloutputs.Freqpll-60.0);
    return true;
  }
  return false;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Epria1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void Epria1::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values+offsetb; // generator array starts from this location

  if(p_mode == XVECTOBUS) {
    p_iabc[0]  = x[0];
    p_iabc[1]  = x[1];
    p_iabc[2]  = x[2];
    p_iLabc[0] = x[3];
    p_iLabc[1] = x[4];
    p_iLabc[2] = x[5];
  } else if(p_mode == XDOTVECTOBUS) {
    p_iLdot[0]  = x[3];
    p_iLdot[1]  = x[4];
    p_iLdot[2]  = x[5];
  } 
}

/**
   Prestep function
*/
void Epria1::preStep(double time ,double timestep)
{
  if(integrationtype != EXPLICIT) {
    return;
  }

  p_vabc[0] = p_va;
  p_vabc[1] = p_vb;
  p_vabc[2] = p_vc;

  double i[3];
  i[0] = p_iabc[0]*sbase/mbase;
  i[1] = p_iabc[1]*sbase/mbase;
  i[2] = p_iabc[2]*sbase/mbase;
  abc2dq0(p_vabc,time,p_Va0,p_vdq0);
  abc2dq0(i,time,p_Va0,p_idq0);

  modelinputs.Va = p_va;
  modelinputs.Vb = p_vb;
  modelinputs.Vc = p_vc;
  
  modelinputs.Ia = p_iabc[0]*sbase/mbase;
  modelinputs.Ib = p_iabc[1]*sbase/mbase;
  modelinputs.Ic = p_iabc[2]*sbase/mbase;

  modelinputs.IaL1 = p_iLabc[0];
  modelinputs.IbL1 = p_iLabc[1];
  modelinputs.IcL1 = p_iLabc[2];

  modelinputs.Time     = time;
  modelinputs.Timestep = timestep;

  int ok = Model_Outputs(&model);

  // Get outputs from the model
  Ed = modeloutputs.Ed;
  Eq = modeloutputs.Eq;
  phi_PLL = modeloutputs.phi_PLL;

}

/**
   Poststep function
*/
void Epria1::postStep(double time)
{
}



/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Epria1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // generator array starts from this location

  if(p_mode == RESIDUAL_EVAL) {
    double e[3];
    double Edq0[3];
    Edq0[0] = Ed;
    Edq0[1] = Eq;
    Edq0[2] = 0.0;

    printf("Time = %lf, Ed = %lf, Eq = %lf, phi = %lf\n",p_time,Ed,Eq,phi_PLL);
    
    dq02abc(Edq0,p_time,phi_PLL,e);

    p_vabc[0] = p_va;
    p_vabc[1] = p_vb;
    p_vabc[2] = p_vc;

    double dvdt[3];
    double Rdamp = modelparams.Rdamp;
    double Cfilt = modelparams.Cfilt/OMEGA_S;

    bus->getVoltageDerivatives(&dvdt[0],&dvdt[1],&dvdt[2]);

    f[0] = (p_iLabc[0] + (zero_Rdamp?0.0:p_va/Rdamp) - Cfilt*dvdt[0])*mbase/sbase - p_iabc[0];
    f[1] = (p_iLabc[1] + (zero_Rdamp?0.0:p_vb/Rdamp) - Cfilt*dvdt[1])*mbase/sbase - p_iabc[1];
    f[2] = (p_iLabc[2] + (zero_Rdamp?0.0:p_vc/Rdamp) - Cfilt*dvdt[2])*mbase/sbase - p_iabc[2];
    
    f[3] = (e[0] - p_Rs*p_iLabc[0] - p_va)/p_L - p_iLdot[0];
    f[4] = (e[1] - p_Rs*p_iLabc[1] - p_vb)/p_L - p_iLdot[1];
    f[5] = (e[2] - p_Rs*p_iLabc[2] - p_vc)/p_L - p_iLdot[2];

  }

}

  /**
   * Return the generator current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
void Epria1::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = p_iabc[0];
  *ib = p_iabc[1];
  *ic = p_iabc[2];
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Epria1::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc;
}


/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values

 Number of non-zero values = 15
 */
int Epria1::matrixNumValues()
{
  int numVals;
  if(integrationtype == IMPLICIT) {
    numVals = 15;
  } else {
    numVals = 15;
  }
  return numVals;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Epria1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  int ia_idx = p_gloc;
  int ib_idx = p_gloc+1;
  int ic_idx = p_gloc+2;

  int iLa_idx = p_gloc+3;
  int iLb_idx = p_gloc+4;
  int iLc_idx = p_gloc+5;

  int va_idx = p_glocvoltage;
  int vb_idx = p_glocvoltage+1;
  int vc_idx = p_glocvoltage+2;

  double Rdamp = modelparams.Rdamp;
  double Cfilt = modelparams.Cfilt/OMEGA_S;

  rows[ctr]   = ia_idx;
  cols[ctr]   = iLa_idx;
  values[ctr]   = 1.0*mbase/sbase;

  rows[ctr+1]   = ia_idx;
  cols[ctr+1]   = va_idx;
  values[ctr+1]   = ((zero_Rdamp?0.0:1/Rdamp) - Cfilt*shift)*mbase/sbase;

  rows[ctr+2]   = ia_idx;
  cols[ctr+2]   = ia_idx;
  values[ctr+2]   = -1.0;

  ctr += 3;

  rows[ctr]   = ib_idx;
  cols[ctr]   = iLb_idx;
  values[ctr]   = 1.0*mbase/sbase;

  rows[ctr+1]   = ib_idx;
  cols[ctr+1]   = vb_idx;
  values[ctr+1]   = ((zero_Rdamp?0.0:1/Rdamp) - Cfilt*shift)*mbase/sbase;

  rows[ctr+2]   = ib_idx;
  cols[ctr+2]   = ib_idx;
  values[ctr+2]   = -1.0;

  ctr += 3;

  rows[ctr]   = ic_idx;
  cols[ctr]   = iLc_idx;
  values[ctr]   = 1.0*mbase/sbase;

  rows[ctr+1]   = ic_idx;
  cols[ctr+1]   = vc_idx;
  values[ctr+1]   = ((zero_Rdamp?0.0:1/Rdamp) - Cfilt*shift)*mbase/sbase;

  rows[ctr+2]   = ic_idx;
  cols[ctr+2]   = ic_idx;
  values[ctr+2]   = -1.0;

  ctr += 3;
  
  rows[ctr]   = iLa_idx;
  cols[ctr]   = iLa_idx;
  values[ctr] = -p_Rs/p_L - shift;

  rows[ctr+1]   = iLa_idx;
  cols[ctr+1]   = va_idx;
  values[ctr+1] = -1.0/p_L;

  ctr += 2;

  rows[ctr]   = iLb_idx;
  cols[ctr]   = iLb_idx;
  values[ctr] = -p_Rs/p_L - shift;

  rows[ctr+1]   = iLb_idx;
  cols[ctr+1]   = vb_idx;
  values[ctr+1] = -1.0/p_L;

  ctr += 2;

  rows[ctr]   = iLc_idx;
  cols[ctr]   = iLc_idx;
  values[ctr] = -p_Rs/p_L - shift;

  rows[ctr+1]   = iLc_idx;
  cols[ctr+1]   = vc_idx;
  values[ctr+1] = -1.0/p_L;

  ctr += 2;

  *nvals = ctr;
}
