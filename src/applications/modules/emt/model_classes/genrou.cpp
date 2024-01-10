#include <genrou.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

/**
   Sauer and Pai - Power System Dynamics and Stability book
   Eqs. 3.148 - 3.159 - Page 42 with the following modifications
   - Stator flux differential equations ignored
   - Additional algebraic equations for machine three-phase currents
   - Speed sensitivity for mechanical torque Tm
*/
Genrou::Genrou(void)
{
  nxgen   = 12; // Number of variables for this model
}

Genrou::~Genrou(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseGeneratorModel
 */
void Genrou::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGenModel::load(data,idx); // load parameters in base generator model

  gridpack::ComplexType Zsource;
  
  // load parameters for the model type
  data->getValue(BUS_NUMBER, &bid);
    // load parameters for the model type
  if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H, &H, idx)) H = 0.0; // H
  if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0, &D, idx)) D = 0.0; // D
  data->getValue(GENERATOR_ZSOURCE,&Zsource,idx);
  Ra = real(Zsource);
  if (!data->getValue(GENERATOR_XD, &Xd, idx)) Xd=0.0; // Xd
  if (!data->getValue(GENERATOR_XQ, &Xq, idx)) Xq=0.0; // Xq
  if (!data->getValue(GENERATOR_XDP, &Xdp, idx)) Xdp=0.0; // Xdp
  if (!data->getValue(GENERATOR_XDPP, &Xdpp, idx)) Xdpp=0.0; // Xdpp
  if (!data->getValue(GENERATOR_XL, &Xl, idx)) Xl=0.0; // Xl
  if (!data->getValue(GENERATOR_TDOP, &Tdop, idx)) Tdop=0.0; // Tdop
  if (!data->getValue(GENERATOR_TDOPP, &Tdopp, idx)) Tdopp=0.0; // Tdopp
  if (!data->getValue(GENERATOR_TQOPP, &Tqopp, idx)) Tqopp=0.0; // Tqopp
  if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.067; // S10 TBD: check parser
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.579; // S12 TBD: check parser
  if (!data->getValue(GENERATOR_XQP, &Xqp, idx)) Xqp=0.0; // Xqp
  if (!data->getValue(GENERATOR_XDPP, &Xqpp, idx)) Xqpp=Xdpp; // Xqpp 
  if (!data->getValue(GENERATOR_TQOP, &Tqop, idx)) Tqop=0.0; // Tqop

  L = Xdpp/OMEGA_S;

}

/**
 * Initialize generator model before calculation
 * @param [output] xin - array where initialized generator variables should be set
 */
void Genrou::init(gridpack::ComplexType* xin)
{
  double IGD,IGQ; // Machine currents in cartesian coordinates
  double Pg, Qg;  // Generator real and reactive power
  double dw0=0.0;  // Initial machine speed deviation
  gridpack::ComplexType *x = xin+offsetb; // generator array starts from this location

  Pg = pg/mbase;
  Qg = qg/mbase;

  VD = p_Vm0*cos(p_Va0);
  VQ = p_Vm0*sin(p_Va0);

  gridpack::ComplexType V = gridpack::ComplexType(VD,VQ);
  gridpack::ComplexType S = gridpack::ComplexType(Pg,Qg);
  gridpack::ComplexType I;
  gridpack::ComplexType Z = gridpack::ComplexType(Ra,Xdpp);
  gridpack::ComplexType Z1 = gridpack::ComplexType(Ra,Xq); // used in delta calculation
  gridpack::ComplexType E,E1;

  // Machine current
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
  
  // Calculation of machine angle delta and coordinate transformation
  // angle theta
  E1 = V + I*Z1;
  delta = arg(E1);
  // theta is behind delta by 90 degrees
  double theta = delta - PI/2.0;

  // Generator internal voltage on network reference frame
  E = V + I*Z;
  double Em = abs(E);
  double Eang = arg(E);

  // Network to machine reference frame transformation
  abc2dq0(vabc,p_time,theta,vdq0);
  abc2dq0(iabc,p_time,theta,idq0);
  
  double Vd, Vq, V0, Id, Iq, I0;

  Vd = vdq0[0];
  Vq = vdq0[1];
  V0 = vdq0[2];

  Id = idq0[0];
  Iq = idq0[1];
  I0 = idq0[2];

  psid = Ra*Iq + Vq;
  psiq = -Ra*Id - Vd;
  psi0 = 0.0;

  dw = dw0;

  psi1d = psid + Xl*Id;
  Eqp   = psi1d + (Xdp - Xl)*Id;

  psi2q = psiq + Xl*Iq;
  Edp = -psi2q - (Xqp - Xl)*Iq;

  TM = psid*Iq - psiq*Id;

  double param, LadIfd;
  double dpsi1ddt;

  param = (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));
  
  dpsi1ddt = psi1d + (Xdp - Xl)*Id - Eqp;

  LadIfd = -Eqp - (Xd - Xdp)*(Id - param*dpsi1ddt);

  Efd = -LadIfd;
  
  // Initialized state variables
  x[0] = psid;
  x[1] = psiq;
  x[2] = psi0;
  x[3] = Eqp;
  x[4] = psi1d;
  x[5] = Edp;
  x[6] = psi2q;
  x[7] = delta;
  x[8] = dw;
  x[9] = iabc[0];
  x[10] = iabc[1];
  x[11] = iabc[2];

}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Genrou::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Genrou::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void Genrou::setValues(gridpack::ComplexType *values)
{
  gridpack::ComplexType *x = values+offsetb; // generator array starts from this location

  if(p_mode == XVECTOBUS) {
    psid  = real(x[0]);
    psiq  = real(x[1]);
    psi0  = real(x[2]);
    Eqp   = real(x[3]);
    psi1d   = real(x[4]);
    Edp     = real(x[5]);
    psi2q = real(x[6]);
    delta = real(x[7]);
    dw = real(x[8]);
    iabc[0] = real(x[9]);
    iabc[1] = real(x[10]);
    iabc[2] = real(x[11]);
  } else if(p_mode == XDOTVECTOBUS) {
    dpsid  = real(x[0]);
    dpsiq  = real(x[1]);
    dpsi0  = real(x[2]);
    dEqp   = real(x[3]);
    dpsi1d   = real(x[4]);
    dEdp     = real(x[5]);
    dpsi2q = real(x[6]);
    ddelta = real(x[7]);
    ddw = real(x[8]);
    diabc[0] = real(x[9]);
    diabc[1] = real(x[10]);
    diabc[2] = real(x[11]);
  } 
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Genrou::vectorGetValues(gridpack::ComplexType *values)
{
  gridpack::ComplexType *f = values+offsetb; // generator array starts from this location

  if(p_mode == RESIDUAL_EVAL) {
    double theta = delta - PI/2.0;

    double tempd1,tempd2,tempq1,tempq2;
    tempd1 = (Xdpp - Xl)/(Xdp - Xl);
    tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
    tempq1 = (Xdpp - Xl)/(Xqp - Xl);
    tempq2 = (Xqp - Xdpp)/(Xqp - Xl);

    vabc[0] = p_va;
    vabc[1] = p_vb;
    vabc[2] = p_vc;

    // Network to machine reference frame transformation
    abc2dq0(vabc,p_time,theta,vdq0);
  
    double Vd, Vq, V0, Id, Iq, I0;

    Vd = vdq0[0];
    Vq = vdq0[1];
    V0 = vdq0[2];

    
    Id = (psid - tempd1*Eqp - tempd2*psi1d)/-Xdpp;
    Iq = (psiq + tempq1*Edp - tempq2*psi2q)/-Xdpp;
    I0 = psi0/-Xl;

    idq0[0] = Id;
    idq0[1] = Iq;
    idq0[2] = I0;

    f[0] = OMEGA_S*(Ra*Id + (1 + dw)*psiq + Vd) - dpsid;
    f[1] = OMEGA_S*(Ra*Iq - (1 + dw)*psid + Vq) - dpsiq;
    f[2] = OMEGA_S*(Ra*I0 + V0) - dpsi0;

    if(hasExciter()) {
      Efd = getExciter()->getFieldVoltage();
    }
		    
    double dpsi1ddt;
    double param1 = (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));

    dpsi1ddt = -psi1d + Eqp - (Xdp - Xl)*Id;

    f[3] = (-Eqp - (Xd - Xdp)*(Id - param1*-dpsi1ddt) + Efd)/Tdop - dEqp;

    f[4] = dpsi1ddt/Tdopp - dpsi1d;

    double dpsi2qdt;
    double param2 = (Xqp - Xdpp)/((Xqp - Xl)*(Xqp - Xl));

    dpsi2qdt = -psi2q - Edp - (Xqp - Xl)*Iq;

    f[5] = (-Edp + (Xq - Xqp)*(Iq - param2*-dpsi2qdt))/Tqop - dEdp;

    f[6] = dpsi2qdt/Tqopp - dpsi2q;

    f[7] = OMEGA_S*dw - ddelta;

    if(hasGovernor()) {
      TM = getGovernor()->getMechanicalPower();
    }
    f[8] = 1 / (2 *H) * ((TM - D*dw)/(1+dw) - (psid*Iq - psiq*Id)) - ddw; 

    double igen[3];
    dq02abc(idq0,p_time,theta,igen);

    f[9] = igen[0]*mbase/sbase - iabc[0];
    f[10] = igen[1]*mbase/sbase - iabc[1];
    f[11] = igen[2]*mbase/sbase - iabc[2];
  }

}

  /**
   * Return the generator current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
void Genrou::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = iabc[0];
  *ib = iabc[1];
  *ic = iabc[2];
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Genrou::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc + 9;
}


/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values

 Non-zero pattern of the Jacobian (x denotes non-zero value)

          psid    psiq    psi0    Eqp    psi1d    Edp    psi2q    delta    omega    ia    ib    ic    va    vb    vc    Efd    PM
 eq. 0 |   x       x               x       x                       x         x                        x      x     x 
 eq. 1 |   x       x                               x       x       x         x                        x      x     x
 eq. 2 |                   x                                       x                                  x      x     x 
 eq. 3 |   x                       x       x                                                                             x                 
 eq. 4 |   x                       x       x         
 eq. 5 |           x                               x       x
 eq. 6 |           x                               x       x
 eq. 7 |                                                           x         x
 eq. 8 |   x       x               x       x       x       x                 x                                                x
 eq. 9 |   x       x       x       x       x       x       x       x                x
 eq.10 |   x       x       x       x       x       x       x       x                      x
 eq.11 |   x       x       x       x       x       x       x       x                            x
 
 Number of non-zero values = 9 + 9 + 5 + 4 + 3 + 3 + 3 + 2 + 7 + 9 + 9 + 9 = 69
 */
int Genrou::matrixNumValues()
{
  int numVals = 69;
  if(hasExciter()) numVals += 1;
  if(hasGovernor()) numVals += 1;

  return numVals;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Genrou::matrixGetValues(int *nvals, gridpack::ComplexType *values, int *rows, int *cols)
{
  int ctr = 0;
  // Set up some indices
  int psid_idx = p_gloc;
  int psiq_idx = p_gloc+1;
  int psi0_idx = p_gloc+2;
  int Eqp_idx  = p_gloc+3;
  int psi1d_idx = p_gloc+4;
  int Edp_idx  = p_gloc+5;
  int psi2q_idx = p_gloc+6;
  int delta_idx = p_gloc+7;
  int dw_idx  = p_gloc+8;
  int ia_idx  = p_gloc+9;
  int ib_idx = p_gloc+10;
  int ic_idx = p_gloc+11;

  double tempd1,tempd2,tempq1,tempq2;
  tempd1 = (Xdpp - Xl)/(Xdp - Xl);
  tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
  tempq1 = (Xdpp - Xl)/(Xqp - Xl);
  tempq2 = (Xqp - Xdpp)/(Xqp - Xl);

  double dId_dpsid, dId_dEqp, dId_dpsi1d;
  double dIq_dpsiq, dIq_dEdp, dIq_dpsi2q;
  double dI0_dpsi0;

  dId_dpsid  = -1/Xdpp;
  dId_dEqp   =  tempd1/Xdpp;
  dId_dpsi1d =  tempd2/Xdpp;

  dIq_dpsiq  = -1/Xdpp;
  dIq_dEdp   = -tempq1/Xdpp;
  dIq_dpsi2q =  tempq2/Xdpp;

  dI0_dpsi0  = -1/Xl;
  
  double vabc[3];
  vabc[0] = p_va;
  vabc[1] = p_vb;
  vabc[2] = p_vc;

  double dvdq0_ddelta[3];
  double Tdq0[3][3],dTdq0_ddelta[3][3];

  double theta = delta - PI/2.0;
  getTdq0(p_time,theta,Tdq0);
  getdTdq0dtheta(p_time,theta,dTdq0_ddelta);

  matvecmult3x3(dTdq0_ddelta,vabc,dvdq0_ddelta);

  // Derivatives w.r.t dpisd_dt
  // psid
  rows[ctr] = psid_idx; cols[ctr] = psid_idx;
  values[ctr] = OMEGA_S*Ra*dId_dpsid - shift;
  
  rows[ctr+1] = psid_idx; cols[ctr+1] = psiq_idx;
  values[ctr+1] = OMEGA_S*(1 + dw);

  ctr += 2;

  rows[ctr] = psid_idx; cols[ctr] = Eqp_idx;
  values[ctr] = OMEGA_S*Ra*dId_dEqp;

  rows[ctr+1] = psid_idx; cols[ctr+1] = psi1d_idx;
  values[ctr+1] = OMEGA_S*Ra*dId_dpsi1d;

  ctr += 2;

  rows[ctr] = psid_idx; cols[ctr] = delta_idx;
  values[ctr] = OMEGA_S*dvdq0_ddelta[0];

  rows[ctr+1] = psid_idx; cols[ctr+1] = dw_idx;
  values[ctr+1] = OMEGA_S*psiq;
  
  ctr += 2;
  
  rows[ctr]   = psid_idx;
  cols[ctr]   = p_glocvoltage;
  values[ctr] = OMEGA_S*Tdq0[0][0];

  rows[ctr+1]   = psid_idx;
  cols[ctr+1]   = p_glocvoltage+1;
  values[ctr+1] = OMEGA_S*Tdq0[0][1];

  rows[ctr+2]   = psid_idx;
  cols[ctr+2]   = p_glocvoltage+2;
  values[ctr+2] = OMEGA_S*Tdq0[0][2];

  ctr += 3;

  assert(ctr == 9);

  // Derivatives w.r.t dpisq_dt

  rows[ctr] = psiq_idx; cols[ctr] = psid_idx;
  values[ctr] = -OMEGA_S*(1 + dw);
  
  rows[ctr+1] = psiq_idx; cols[ctr+1] = psiq_idx;
  values[ctr+1] = OMEGA_S*Ra*dIq_dpsiq - shift;

  ctr += 2;

  rows[ctr] = psiq_idx; cols[ctr] = Edp_idx;
  values[ctr] = OMEGA_S*Ra*dIq_dEdp;

  rows[ctr+1] = psiq_idx; cols[ctr+1] = psi2q_idx;
  values[ctr+1] = OMEGA_S*Ra*dIq_dpsi2q;

  ctr += 2;

  rows[ctr] = psiq_idx; cols[ctr] = delta_idx;
  values[ctr] = OMEGA_S*dvdq0_ddelta[1];

  rows[ctr+1] = psiq_idx; cols[ctr+1] = dw_idx;
  values[ctr+1] = -OMEGA_S*psid;
  
  ctr += 2;
  
  rows[ctr]   = psiq_idx;
  cols[ctr]   = p_glocvoltage;
  values[ctr] = OMEGA_S*Tdq0[1][0];

  rows[ctr+1]   = psiq_idx;
  cols[ctr+1]   = p_glocvoltage+1;
  values[ctr+1] = OMEGA_S*Tdq0[1][1];

  rows[ctr+2]   = psiq_idx;
  cols[ctr+2]   = p_glocvoltage+2;
  values[ctr+2] = OMEGA_S*Tdq0[1][2];

  ctr += 3;

  assert(ctr == 18);

  // Derivatives w.r.t dpsi0_dt

  rows[ctr] = psi0_idx;  cols[ctr] = psi0_idx;
  values[ctr] = OMEGA_S*Ra*dI0_dpsi0 - shift;

  ctr += 1;

  rows[ctr] = psi0_idx; cols[ctr] = delta_idx;
  values[ctr] = OMEGA_S*dvdq0_ddelta[2];

  ctr += 1;

  rows[ctr]   = psi0_idx;
  cols[ctr]   = p_glocvoltage;
  values[ctr] = OMEGA_S*Tdq0[2][0];

  rows[ctr+1]   = psi0_idx;
  cols[ctr+1]   = p_glocvoltage+1;
  values[ctr+1] = OMEGA_S*Tdq0[2][1];

  rows[ctr+2]   = psi0_idx;
  cols[ctr+2]   = p_glocvoltage+2;
  values[ctr+2] = OMEGA_S*Tdq0[2][2];

  ctr += 3;
  
  assert(ctr == 23);

  // Derivatives w.r.t. dEqp_dt
  double dpsi1ddt_dpsi1d, dpsi1ddt_dEqp, dpsi1ddt_dpsid;
  double param1 = (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));

  dpsi1ddt_dpsid = -(Xdp - Xl)*dId_dpsid;
  dpsi1ddt_dEqp  = 1.0 -(Xdp - Xl)*dId_dEqp;
  dpsi1ddt_dpsi1d = -1.0 -(Xdp - Xl)*dId_dpsi1d;
  
  rows[ctr] = Eqp_idx;  cols[ctr] = psid_idx;
  values[ctr] = (-(Xd - Xdp)*(dId_dpsid - param1*-dpsi1ddt_dpsid))/Tdop;

  rows[ctr+1] = Eqp_idx; cols[ctr+1] = Eqp_idx;
  values[ctr+1] = (-1.0 -(Xd - Xdp)*(dId_dEqp - param1*(-dpsi1ddt_dEqp)))/Tdop -shift;

  rows[ctr+2] = Eqp_idx; cols[ctr+2] = psi1d_idx;
  values[ctr+2] = (-(Xd - Xdp)*(dId_dpsi1d -param1*-dpsi1ddt_dpsi1d))/Tdop;

  ctr += 3;
  if(hasExciter()) {
    int Efd_idx;
    double Efd;
    Efd = getExciter()->getFieldVoltage(&Efd_idx);
    rows[ctr] = Eqp_idx; cols[ctr] = Efd_idx;
    values[ctr] = 1.0/Tdop;
    ctr += 1;
  }

  // Derivative of dpsi1d_dt
  rows[ctr] = psi1d_idx;  cols[ctr] = psid_idx;
  values[ctr] = (dpsi1ddt_dpsid)/Tdopp;

  rows[ctr+1] = psi1d_idx; cols[ctr+1] = Eqp_idx;
  values[ctr+1] = (dpsi1ddt_dEqp)/Tdopp;

  rows[ctr+2] = psi1d_idx; cols[ctr+2] = psi1d_idx;
  values[ctr+2] = (dpsi1ddt_dpsi1d)/Tdopp - shift;

  ctr += 3;

  // Derivatives w.r.t. dEdp_dt
  double dpsi2qdt_dpsi2q, dpsi2qdt_dEdp, dpsi2qdt_dpsiq;
  double param2 = (Xqp - Xdpp)/((Xqp - Xl)*(Xqp - Xl));

  dpsi2qdt_dpsiq = -(Xqp - Xl)*dIq_dpsiq;
  dpsi2qdt_dEdp  = -1.0 -(Xqp - Xl)*dIq_dEdp;
  dpsi2qdt_dpsi2q = -1.0 -(Xqp - Xl)*dIq_dpsi2q;
  
  rows[ctr] = Edp_idx;  cols[ctr] = psiq_idx;
  values[ctr] = ((Xq - Xqp)*(dIq_dpsiq - param2*-dpsi2qdt_dpsiq))/Tqop;

  rows[ctr+1] = Edp_idx; cols[ctr+1] = Edp_idx;
  values[ctr+1] = (-1.0 +(Xq - Xqp)*(dIq_dEdp - param2*(-dpsi2qdt_dEdp)))/Tqop - shift;

  rows[ctr+2] = Edp_idx; cols[ctr+2] = psi2q_idx;
  values[ctr+2] = ((Xq - Xqp)*(dIq_dpsi2q -param2*-dpsi2qdt_dpsi2q))/Tqop;

  ctr += 3;

  // Derivative of dpsi2q_dt
  rows[ctr] = psi2q_idx;  cols[ctr] = psiq_idx;
  values[ctr] = (dpsi2qdt_dpsiq)/Tqopp;

  rows[ctr+1] = psi2q_idx; cols[ctr+1] = Edp_idx;
  values[ctr+1] = (dpsi2qdt_dEdp)/Tqopp;

  rows[ctr+2] = psi2q_idx; cols[ctr+2] = psi2q_idx;
  values[ctr+2] = (dpsi2qdt_dpsi2q)/Tqopp - shift;

  ctr += 3;

  // Derivative of ddelta_dt
  rows[ctr] = delta_idx; cols[ctr] = delta_idx;
  values[ctr] = -shift;

  rows[ctr+1] = delta_idx; cols[ctr+1] = dw_idx;
  values[ctr+1] = OMEGA_S;

  ctr += 2;

  // derivative of ddw_dt
  double Minv = 1 / (2*H);
  double Id,Iq,I0;

  Id = (psid - tempd1*Eqp - tempd2*psi1d)/-Xdpp;
  Iq = (psiq + tempq1*Edp - tempq2*psi2q)/-Xdpp;
  I0 = psi0/-Xl;

  rows[ctr] = dw_idx; cols[ctr] = psid_idx;
  values[ctr] = Minv*-(Iq - psiq*dId_dpsid);

  rows[ctr+1] = dw_idx; cols[ctr+1] = psiq_idx;
  values[ctr+1] = Minv*-(psid*dIq_dpsiq - Id);

  rows[ctr+2] = dw_idx; cols[ctr+2] = Eqp_idx;
  values[ctr+2] = Minv*-(-psiq*dId_dEqp);

  rows[ctr+3] = dw_idx; cols[ctr+3] = psi1d_idx;
  values[ctr+3] = Minv*-(-psiq*dId_dpsi1d);

  rows[ctr+4] = dw_idx; cols[ctr+4] = Edp_idx;
  values[ctr+4] = Minv*-(psid*dIq_dEdp);

  rows[ctr+5] = dw_idx; cols[ctr+5] = psi2q_idx;
  values[ctr+5] = Minv*-(psid*dIq_dpsi2q);

  rows[ctr+6] = dw_idx; cols[ctr+6] = dw_idx;
  values[ctr+6] = -D*(1/(1+dw) - dw/((1+dw)*(1+dw))) -shift;

  ctr += 7;
  if(hasGovernor()) {
    // Partial derivatives w.r.t Governor
  }
  
  // derivative of currents iabc

  rows[ctr]   = ia_idx; cols[ctr]   = psid_idx;
  rows[ctr+1] = ia_idx; cols[ctr+1] = psiq_idx;
  rows[ctr+2] = ia_idx; cols[ctr+2] = psi0_idx;
  rows[ctr+3] = ia_idx; cols[ctr+3] = Eqp_idx;
  rows[ctr+4] = ia_idx; cols[ctr+4] = psi1d_idx;
  rows[ctr+5] = ia_idx; cols[ctr+5] = Edp_idx;
  rows[ctr+6] = ia_idx; cols[ctr+6] = psi2q_idx;
  rows[ctr+7] = ia_idx; cols[ctr+7] = delta_idx;
  rows[ctr+8] = ia_idx; cols[ctr+8] = ia_idx;

  double Tdq0inv[3][3],dTdq0inv_ddelta[3][3];
  double scal = sbase/mbase;

  getTdq0inv(p_time,theta,Tdq0inv);
  getdTdq0invdtheta(p_time,theta,dTdq0inv_ddelta);

  values[ctr]   = scal*(Tdq0inv[0][0]*dId_dpsid);
  values[ctr+1] = scal*(Tdq0inv[0][1]*dIq_dpsiq);
  values[ctr+2] = scal*(Tdq0inv[0][2]*dI0_dpsi0);

  values[ctr+3] = scal*(Tdq0inv[0][0]*dId_dEqp);
  values[ctr+4] = scal*(Tdq0inv[0][0]*dId_dpsi1d);

  values[ctr+5] = scal*(Tdq0inv[0][1]*dIq_dEdp);
  values[ctr+6] = scal*(Tdq0inv[0][1]*dIq_dpsi2q);
  
  values[ctr+7] = scal*(dTdq0inv_ddelta[0][0]*Id + dTdq0inv_ddelta[0][1]*Iq + dTdq0inv_ddelta[0][2]*I0);
			
  values[ctr+8] = -1.0;
  
  ctr += 9;

  rows[ctr]   = ib_idx; cols[ctr]   = psid_idx;
  rows[ctr+1] = ib_idx; cols[ctr+1] = psiq_idx;
  rows[ctr+2] = ib_idx; cols[ctr+2] = psi0_idx;
  rows[ctr+3] = ib_idx; cols[ctr+3] = Eqp_idx;
  rows[ctr+4] = ib_idx; cols[ctr+4] = psi1d_idx;
  rows[ctr+5] = ib_idx; cols[ctr+5] = Edp_idx;
  rows[ctr+6] = ib_idx; cols[ctr+6] = psi2q_idx;
  rows[ctr+7] = ib_idx; cols[ctr+7] = delta_idx;
  rows[ctr+8] = ib_idx; cols[ctr+8] = ib_idx;

  values[ctr]   = scal*(Tdq0inv[1][0]*dId_dpsid);
  values[ctr+1] = scal*(Tdq0inv[1][1]*dIq_dpsiq);
  values[ctr+2] = scal*(Tdq0inv[1][2]*dI0_dpsi0);

  values[ctr+3] = scal*(Tdq0inv[1][0]*dId_dEqp);
  values[ctr+4] = scal*(Tdq0inv[1][0]*dId_dpsi1d);

  values[ctr+5] = scal*(Tdq0inv[1][1]*dIq_dEdp);
  values[ctr+6] = scal*(Tdq0inv[1][1]*dIq_dpsi2q);
  
  values[ctr+7] = scal*(dTdq0inv_ddelta[1][0]*Id + dTdq0inv_ddelta[1][1]*Iq + dTdq0inv_ddelta[1][2]*I0);
			
  values[ctr+8] = -1.0;
  
  ctr += 9;

  rows[ctr]   = ic_idx; cols[ctr]   = psid_idx;
  rows[ctr+1] = ic_idx; cols[ctr+1] = psiq_idx;
  rows[ctr+2] = ic_idx; cols[ctr+2] = psi0_idx;
  rows[ctr+3] = ic_idx; cols[ctr+3] = Eqp_idx;
  rows[ctr+4] = ic_idx; cols[ctr+4] = psi1d_idx;
  rows[ctr+5] = ic_idx; cols[ctr+5] = Edp_idx;
  rows[ctr+6] = ic_idx; cols[ctr+6] = psi2q_idx;
  rows[ctr+7] = ic_idx; cols[ctr+7] = delta_idx;
  rows[ctr+8] = ic_idx; cols[ctr+8] = ic_idx;

  values[ctr]   = scal*(Tdq0inv[2][0]*dId_dpsid);
  values[ctr+1] = scal*(Tdq0inv[2][1]*dIq_dpsiq);
  values[ctr+2] = scal*(Tdq0inv[2][2]*dI0_dpsi0);

  values[ctr+3] = scal*(Tdq0inv[2][0]*dId_dEqp);
  values[ctr+4] = scal*(Tdq0inv[2][0]*dId_dpsi1d);

  values[ctr+5] = scal*(Tdq0inv[2][1]*dIq_dEdp);
  values[ctr+6] = scal*(Tdq0inv[2][1]*dIq_dpsi2q);
  
  values[ctr+7] = scal*(dTdq0inv_ddelta[2][0]*Id + dTdq0inv_ddelta[2][1]*Iq + dTdq0inv_ddelta[2][2]*I0);
			
  values[ctr+8] = -1.0;
  
  ctr += 9;

  *nvals = ctr;
}


/**
 * Set Jacobian values
 * @param values a 2-d array of Jacobian block for the bus
 */
bool Genrou::setJacobian(gridpack::ComplexType **values)
{
  return true;
}

/**
 * Returns the initial field voltage (Efd(t0))
 * @param [out] Efd0 - Initial field voltage
 */
double Genrou::getInitialFieldVoltage()
{
  return Efd;
}
