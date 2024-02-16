#include <gencls.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Gencls::Gencls(void)
{
  p_delta = 0.0;
  p_dw    = 0.0;
  p_deltadot = 0.0;
  p_dwdot    = 0.0;
  p_Rs    = 0.0;
  p_L     = 1.0/OMEGA_S;
  p_H     = 0.0;
  p_D     = 0.0;
  p_Ep    = 0.0;
  p_Pm    = 0.0;
  p_Xdp   = 0.0;
  
  nxgen   = 5; // Number of variables for this model when integration type is implicit or implicit-explicit
}

void Gencls::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgen = 3;
  *nvar = nxgen;
}

Gencls::~Gencls(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseGeneratorModel
 */
void Gencls::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGenModel::load(data,idx); // load parameters in base generator model
  gridpack::ComplexType Zsource;

  // load parameters for the model type
  data->getValue(BUS_NUMBER, &bid);
  data->getValue(GENERATOR_ZSOURCE,&Zsource,idx);
  p_Rs = real(Zsource);
  p_Xdp = imag(Zsource);
  //  data->getValue(GENERATOR_RESISTANCE,&p_Rs,idx);
  //  data->getValue(GENERATOR_TRANSIENT_REACTANCE,&p_Xdp,idx);
  data->getValue(GENERATOR_INERTIA_CONSTANT_H,&p_H,idx);
  data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&p_D,idx);

  p_L = p_Xdp/OMEGA_S;
}

/**
 * Initialize generator model before calculation
 * @param [output] xin - array where initialized generator variables should be set
 */
void Gencls::init(gridpack::RealType* xin)
{
  double IGD,IGQ; // Machine currents in cartesian coordinates
  double Pg, Qg;  // Generator real and reactive power
  double dw=0.0;  // Initial machine speed deviation
  gridpack::RealType *x = xin+offsetb; // generator array starts from this location

  Pg = pg/sbase;
  Qg = qg/sbase;

  VD = p_Vm0*cos(p_Va0);
  VQ = p_Vm0*sin(p_Va0);

  gridpack::ComplexType V = gridpack::ComplexType(VD,VQ);
  gridpack::ComplexType S = gridpack::ComplexType(Pg,Qg);
  gridpack::ComplexType I;
  gridpack::ComplexType Z = gridpack::ComplexType(p_Rs,p_Xdp);
  gridpack::ComplexType E;

  I = conj(S/V);
  E = V + I*Z;

  double delta = arg(E);
  double Im = abs(I);
  double Ia = arg(I);

  p_Ep = abs(E);
  p_Pm = Pg;

  double ia,ib,ic;

  ia = Im*sin(Ia);
  ib = Im*sin(Ia - 2*PI/3.0);
  ic = Im*sin(Ia + 2*PI/3.0);

  p_iabc[0] = ia;
  p_iabc[1] = ib;
  p_iabc[2] = ic;

  double p_Pm = (p_va*ia + p_vb*ib + p_vc*ic)*(2.0/3.0);
	
  x[0] = ia;
  x[1] = ib;
  x[2] = ic;
  
  if(integrationtype != EXPLICIT) {
    x[3] = delta;
    x[4] = dw;
  } else {
    p_delta = delta;
    p_dw    = dw;
  }
}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Gencls::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Gencls::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void Gencls::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values+offsetb; // generator array starts from this location

  if(p_mode == XVECTOBUS) {
    p_iabc[0]  = x[0];
    p_iabc[1]  = x[1];
    p_iabc[2]  = x[2];
    if(integrationtype != EXPLICIT) {
      p_delta    = x[3];
      p_dw       = x[4];
    }
  } else if(p_mode == XDOTVECTOBUS) {
    p_idot[0]  = x[0];
    p_idot[1]  = x[1];
    p_idot[2]  = x[2];
    if(integrationtype != EXPLICIT) {
      p_deltadot = x[3];
      p_dwdot    = x[4];
    }
  } 
}

/**
   Prestep function
*/
void Gencls::preStep(double time ,double timestep)
{
  if(integrationtype != EXPLICIT) {
    return;
  }
  double Pe;
  double ddelta_dt,ddw_dt;

  p_vabc[0] = p_va;
  p_vabc[1] = p_vb;
  p_vabc[2] = p_vc;

  abc2dq0(p_vabc,p_time,p_delta,p_vdq0);
  abc2dq0(p_iabc,p_time,p_delta,p_idq0);

  Pe = p_vdq0[0]*p_idq0[0] + p_vdq0[1]*p_idq0[1] + p_vdq0[2]*p_idq0[2];
    
  ddelta_dt = p_dw/OMEGA_S;
  ddw_dt    = (p_Pm - Pe - p_D*p_dw)/(2*p_H);

  p_delta = p_delta + timestep*ddelta_dt;
  p_dw    = p_dw    + timestep*ddw_dt;

}

/**
   Poststep function
*/
void Gencls::postStep(double time)
{
}



/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Gencls::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // generator array starts from this location

  if(p_mode == RESIDUAL_EVAL) {
    double Pe;
    double e[3];

    e[0] = p_Ep*sin(OMEGA_S*p_time + p_delta);
    e[1] = p_Ep*sin(OMEGA_S*p_time + p_delta - TWOPI_OVER_THREE);
    e[2] = p_Ep*sin(OMEGA_S*p_time + p_delta + TWOPI_OVER_THREE);

    p_vabc[0] = p_va;
    p_vabc[1] = p_vb;
    p_vabc[2] = p_vc;

    abc2dq0(p_vabc,p_time,p_delta,p_vdq0);
    abc2dq0(p_iabc,p_time,p_delta,p_idq0);

    Pe = p_vdq0[0]*p_idq0[0] + p_vdq0[1]*p_idq0[1] + p_vdq0[2]*p_idq0[2];
    
    // Generator equations
    if(abs(p_L) > 1e-6) {
      // f = di_dt - idot => L^-1*(e - R*i - v) - idot
      f[0] = (e[0] - p_Rs*p_iabc[0] - p_va)/p_L - p_idot[0];
      f[1] = (e[1] - p_Rs*p_iabc[1] - p_vb)/p_L - p_idot[1];
      f[2] = (e[2] - p_Rs*p_iabc[2] - p_vc)/p_L - p_idot[2];
    } else {
      f[0] = (e[0] - p_Rs*p_iabc[0] - p_va);
      f[1] = (e[1] - p_Rs*p_iabc[1] - p_vb);
      f[2] = (e[2] - p_Rs*p_iabc[2] - p_vc);
    }

    if(integrationtype != EXPLICIT) {
      f[3] = p_dw/OMEGA_S - p_deltadot;
      f[4]    = (p_Pm - Pe - p_D*p_dw)/(2*p_H) - p_dwdot;
    }
  }

}

  /**
   * Return the generator current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
void Gencls::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = p_iabc[0];
  *ib = p_iabc[1];
  *ic = p_iabc[2];
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Gencls::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc;
}


/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values

 Non-zero pattern of the Jacobian (x denotes non-zero value)

          ia   ib   ic    delta    omega    va    vb    vc
 eq. 1 |   x                x                x
 eq. 2 |        x           x                      x
 eq. 3 |             x      x                            x
 eq. 4 |                    x       x    
 eq. 5 |   x    x    x      x       x        x     x     x

 Number of non-zero values = 19
 */
int Gencls::matrixNumValues()
{
  int numVals;
  if(integrationtype == IMPLICIT) {
    numVals = 19;
  } else {
    numVals = 6;
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
void Gencls::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  int ia_idx = p_gloc;
  int ib_idx = p_gloc+1;
  int ic_idx = p_gloc+2;

  int delta_idx = p_gloc+3;
  int dw_idx = p_gloc+4;

  int va_idx = p_glocvoltage;
  int vb_idx = p_glocvoltage+1;
  int vc_idx = p_glocvoltage+2;

  // partial derivatives w.r.t. f[0]-f[2]
  if(abs(p_L) > 1e-6) {
    // f[0] partial derivative w.r.t. delta
    if(integrationtype == IMPLICIT) {
      rows[ctr]   = ia_idx;
      cols[ctr]   = delta_idx;
      values[ctr] = p_Ep*cos(OMEGA_S*p_time + p_delta)/p_L;

      ctr += 1;
    }
    rows[ctr]   = ia_idx;
    cols[ctr]   = ia_idx;
    values[ctr] = -p_Rs/p_L - shift;

    rows[ctr+1]   = ia_idx;
    cols[ctr+1]   = va_idx;
    values[ctr+1] = -1.0/p_L;

    ctr += 2;

    if(integrationtype == IMPLICIT) {
      // f[1] partial derivative w.r.t. delta
      rows[ctr]   = ib_idx;
      cols[ctr]   = delta_idx;
      values[ctr] = p_Ep*cos(OMEGA_S*p_time + p_delta - TWOPI_OVER_THREE)/p_L;

      ctr += 1;
    }
    
    rows[ctr]   = ib_idx;
    cols[ctr]   = ib_idx;
    values[ctr] = -p_Rs/p_L - shift;

    rows[ctr+1]   = ib_idx;
    cols[ctr+1]   = vb_idx;
    values[ctr+1] = -1.0/p_L;

    ctr += 2;

    // f[2] partial derivative w.r.t. delta
    if(integrationtype == IMPLICIT) {
      rows[ctr]   = ic_idx;
      cols[ctr]   = delta_idx;
      values[ctr] = p_Ep*cos(OMEGA_S*p_time + p_delta + TWOPI_OVER_THREE)/p_L;

      ctr += 1;
    }
    
    rows[ctr]   = ic_idx;
    cols[ctr]   = ic_idx;
    values[ctr] = -p_Rs/p_L - shift;

    rows[ctr+1]   = ic_idx;
    cols[ctr+1]   = vc_idx;
    values[ctr+1] = -1.0/p_L;

    ctr += 2;
  } else {
    if(integrationtype == IMPLICIT) {
      // f[0] partial derivative w.r.t. delta
      rows[ctr]   = ia_idx;
      cols[ctr]   = delta_idx;
      values[ctr] = p_Ep*cos(OMEGA_S*p_time + p_delta);
      ctr += 1;
    }

    rows[ctr]   = ia_idx;
    cols[ctr]   = ia_idx;
    values[ctr] = -p_Rs;

    rows[ctr+1]   = ia_idx;
    cols[ctr+1]   = va_idx;
    values[ctr+1] = -1.0;

    ctr += 2;

    if(integrationtype == IMPLICIT) {
      // f[1] partial derivative w.r.t. delta
      rows[ctr]   = ib_idx;
      cols[ctr]   = delta_idx;
      values[ctr] = p_Ep*cos(OMEGA_S*p_time + p_delta - TWOPI_OVER_THREE);
      ctr += 1;
    }

    rows[ctr]   = ib_idx;
    cols[ctr]   = ib_idx;
    values[ctr] = -p_Rs;

    rows[ctr+1]   = ib_idx;
    cols[ctr+1]   = vb_idx;
    values[ctr+1] = -1.0;

    ctr += 2;

    if(integrationtype == IMPLICIT) {
      // f[2] partial derivative w.r.t. delta
      rows[ctr]   = ic_idx;
      cols[ctr]   = delta_idx;
      values[ctr] = p_Ep*cos(OMEGA_S*p_time + p_delta + TWOPI_OVER_THREE);
      ctr += 1;
    }
    
    rows[ctr]   = ic_idx;
    cols[ctr]   = ic_idx;
    values[ctr] = -p_Rs;

    rows[ctr+1]   = ic_idx;
    cols[ctr+1]   = vc_idx;
    values[ctr+1] = -1.0;

    ctr += 2;
  }

  if(integrationtype == IMPLICIT) {
    // partial derivatives w.r.t f[3]
    rows[ctr]   = delta_idx;
    cols[ctr]   = delta_idx;
    values[ctr] = -shift;
    
    rows[ctr+1]   = delta_idx;
    cols[ctr+1]   = dw_idx;
    values[ctr+1] = 1/OMEGA_S;
    
    ctr += 2;
    
    // partial derivatives w.r.t. f[4]
    // df[1]_ddw
    rows[ctr]   = dw_idx;
    cols[ctr]   = dw_idx;
    values[ctr] = -p_D/(2*p_H) - shift;
    ctr++;
    
    // f[4] = -1/(2*H)*Pe = -1/(2*H)*vdq0^T*idq0
    // df[4]_ddelta = -1/(2*H)*[vdq0^T*didq0_ddelta + (dvdq0_ddelta)^T*idq0]
    //              = -1/(2*H)*[vdq0^T*dTdq0_ddelta*iabc + (dTdq0_ddelta*vabc)^T*idq0]
    double Tdq0[3][3],dTdq0_ddelta[3][3];
    double omegat = OMEGA_S*p_time + p_delta;
    double omegat_minus = omegat - 2.0*PI/3.0;
    double omegat_plus = omegat + 2.0*PI/3.0;
    
    Tdq0[0][0] = TWO_OVER_THREE*sin(omegat); Tdq0[0][1] = TWO_OVER_THREE*sin(omegat_minus); Tdq0[0][2] = TWO_OVER_THREE*sin(omegat_plus);
    Tdq0[1][0] = TWO_OVER_THREE*cos(omegat); Tdq0[1][1] = TWO_OVER_THREE*cos(omegat_minus); Tdq0[1][2] = TWO_OVER_THREE*cos(omegat_plus);
    Tdq0[2][0] = Tdq0[2][1] = Tdq0[2][2] = TWO_OVER_THREE*0.5;
    
    dTdq0_ddelta[0][0] = TWO_OVER_THREE*cos(omegat); dTdq0_ddelta[0][1] = TWO_OVER_THREE*cos(omegat_minus); dTdq0_ddelta[0][2] = TWO_OVER_THREE*cos(omegat_plus);
    dTdq0_ddelta[1][0] = TWO_OVER_THREE*-sin(omegat); dTdq0_ddelta[1][1] = TWO_OVER_THREE*-sin(omegat_minus); dTdq0_ddelta[1][2] = TWO_OVER_THREE*-sin(omegat_plus);
    Tdq0[2][0] = Tdq0[2][1] = Tdq0[2][2] = 0.0;
    
    double dvdq0_ddelta[3];
    double didq0_ddelta[3];
    double dPe_ddelta1,dPe_ddelta2;
    
    matvecmult3x3(dTdq0_ddelta,p_vabc,dvdq0_ddelta);
    matvecmult3x3(dTdq0_ddelta,p_iabc,didq0_ddelta);
    
    vecdot3(p_vdq0,didq0_ddelta,&dPe_ddelta1);
    vecdot3(p_idq0,dvdq0_ddelta,&dPe_ddelta2);
    
    rows[ctr] = dw_idx;
    cols[ctr] = delta_idx;
    values[ctr] = -1.0/(2*p_H)*(dPe_ddelta1 + dPe_ddelta2);
    
    ctr++;
    
    // df[1]_dvabc = -1/(2*H)*[(dvdq0_dvabc)^T*idq0]
    //             = -1/(2*H)*[d(vabc^T*Tdq0^T)*idq0]
    //              = -1/(2*H)*[Tdq0^T*idq0]
    double df4_dvabc[3];
    
    scaledmattransposevecmult3x3(Tdq0,p_idq0,df4_dvabc,-1.0/(2*p_H));
    
    rows[ctr]   = dw_idx;
    cols[ctr]   = va_idx;
    values[ctr] = df4_dvabc[0];
    
    rows[ctr+1]   = dw_idx;
    cols[ctr+1]   = vb_idx;
    values[ctr+1] = df4_dvabc[1];
    
    rows[ctr+2]   = dw_idx;
    cols[ctr+2]   = vc_idx;
    values[ctr+2] = df4_dvabc[2];
    
    ctr += 3;
    
    // df[1]_diabc = -1/(2*H)*[vdq0^T*Tdq0]
    //              = -1/(2*H)*[vdq0^T*Tdq0]]
    double df4_diabc[3];
    
    scaledvec3multmat3x3(p_vdq0,Tdq0,df4_diabc,-1.0/(2*p_H));
    
    rows[ctr]   = dw_idx;
    cols[ctr]   = ia_idx;
    values[ctr] = df4_diabc[0];
    
    rows[ctr+1]   = dw_idx;
    cols[ctr+1]   = ib_idx;
    values[ctr+1] = df4_diabc[1];
    
    rows[ctr+2]   = dw_idx;
    cols[ctr+2]   = ic_idx;
    values[ctr+2] = df4_diabc[2];
    
    ctr += 3;
  }
  
  *nvals = ctr;
}
