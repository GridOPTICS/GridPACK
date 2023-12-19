#include <genrou.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Genrou::Genrou(void)
{
  delta = 0.0;
  dw    = 0.0;
  Eqp   = 0.0;
  Psidp = 0.0;
  Psiqp = 0.0;
  Edp   = 0.0;

  ddelta = 0.0;
  ddw    = 0.0;
  dEqp   = 0.0;
  dPsidp = 0.0;
  dPsiqp = 0.0;
  dEdp   = 0.0;

  nxgen   = 9; // Number of variables for this model

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

  // load parameters for the model type
  data->getValue(BUS_NUMBER, &bid);
    // load parameters for the model type
  if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H, &H, idx)) H = 0.0; // H
  if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0, &D, idx)) D = 0.0; // D
  if (!data->getValue(GENERATOR_RESISTANCE, &Ra, idx)) Ra=0.0; // Ra
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

  // Machine mechanical power input
  Pmech = Pg;

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

  iabc[0] = ia;
  iabc[1] = ib;
  iabc[2] = ic;

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
  
  double Eppr = real(E);
  double Eppi = imag(E);

  double eabc[3]; // Internal voltage in network reference frame
  double edq0[3]; // Internal voltage in dq axis reference frame

  eabc[0] = Em*sin(OMEGA_S*p_time + Eang);
  eabc[1] = Em*sin(OMEGA_S*p_time + Eang - TWOPI_OVER_THREE);
  eabc[2] = Em*sin(OMEGA_S*p_time + Eang + TWOPI_OVER_THREE);

  abc2dq0(eabc,p_time,theta,edq0);

  double Vd, Vq, Id, Iq;

  Vd = vdq0[0];
  Vq = vdq0[1];
  Id = idq0[0];
  Iq = idq0[1];
  
  double Psidpp = Vq;
  double Psiqpp = -Vd;
  double Psipp = sqrt(Psidpp*Psidpp + Psiqpp*Psiqpp);

  double Psid = Psidpp - Id*Xdpp;
  double Psiq = Psiqpp - Iq*Xdpp;

  Psidp = Psidpp - (Xdpp - Xl)*Id;
  Eqp   = Psidp  + (Xdp - Xl)*Id;

  Edp = (Xq - Xqp)*Iq; // - (Xq - Xl)/(Xd - Xl)*Psiqpp*Sat(Psidpp,Psiqpp)/Psipp;
  Psiqp = Edp + (Xqp - Xl)*Iq;

  Efd = Eqp + (Xd - Xdp)*Id; //  + Psidpp*Sat(Psidpp,Psiqpp)/Psipp;;

  LadIfd = Efd;

  double Telec = Psidpp*Iq - Psiqpp*Id;

  // Initialized state variables
  x[0] = delta;
  x[1] = dw;
  x[2] = Eqp;
  x[3] = Psidp;
  x[4] = Psiqp;
  x[5] = Edp;
  x[6] = iabc[0];
  x[7] = iabc[1];
  x[8] = iabc[2];

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
  } else if(p_mode == XDOTVECTOBUS) {
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
  *i_gloc = p_gloc + 6;
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



