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

         delta    omega    ia    ib    ic    va    vb    vc
 eq. 1 |   x       x    
 eq. 2 |   x       x        x     x     x     x     x     x
 eq. 3 |   x                x                 x
 eq. 4 |   x                      x                 x
 eq. 5 |   x                            x                 x

 Number of non-zero values = 19
 */
int Genrou::matrixNumValues()
{
  int numVals = 19;

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

  if(abs(L) > 1e-6) {
  } else {
  }

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
