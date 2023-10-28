#include <constantimpedance.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Constantimpedance::Constantimpedance(void)
{
  nxload   = 3; // Number of variables for this model (assuming the reactive part of load exists)
}

Constantimpedance::~Constantimpedance(void)
{
}

/**
 * Load parameters from DataCollection object into load model
 * @param data collection of load parameters from input files
 * @param index of load on bus
 * TODO: might want to move this functionality to BaseLoadModel
 */
void Constantimpedance::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTLoadModel::load(data,idx); // load parameters in base load model

  pl /= sbase; // per unit conversion
  ql /= sbase; // per unit conversion
  if(ql == 0.0) {
    // No inductor current
    nxload = 0;
  } else {
    nxload = 3;
  }

}

/**
 * Initialize load model before calculation
 * @param [output] xin - array where initialized load variables should be set
 */
void Constantimpedance::init(gridpack::RealType* xin)
{
  double VD,VQ;
  double Yp,Yq;
  double R,L;
  double Im,Ia;
  gridpack::RealType *x = xin + offsetb;

  VD = p_Vm0*cos(p_Va0);
  VQ = p_Vm0*sin(p_Va0);

  gridpack::ComplexType V = gridpack::ComplexType(VD,VQ);
  gridpack::ComplexType S = gridpack::ComplexType(pl,ql);
  gridpack::ComplexType I = conj(S/V);
  gridpack::ComplexType Z = V/I;
  gridpack::ComplexType Y = 1.0/Z;

  Im = abs(I);
  Ia = atan2(imag(I),real(I));

  R = real(Z);
  L = imag(Z)/OMEGA_S;

  p_R[0] = p_R[1] = p_R[2] = R;
  p_L[0] = p_L[1] = p_L[2] = L;

  x[0] = p_i[0] = Im*sin(Ia);
  x[1] = p_i[1] = Im*sin(Ia-2*PI/3.0);
  x[2] = p_i[2] = Im*sin(Ia+2*PI/3.0);

}

/**
 * Write output from loads to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Constantimpedance::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out load state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Constantimpedance::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto loads
 * @param values array containing load state variables
*/
void Constantimpedance::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values+offsetb; // load array starts from this location

  if(p_mode == XVECTOBUS) {
    p_i[0]  = x[0];
    p_i[1]  = x[1];
    p_i[2]  = x[2];
  } else if(p_mode == XDOTVECTOBUS) {
    p_idot[0]  = x[0];
    p_idot[1]  = x[1];
    p_idot[2]  = x[2];
  } 

}

/**
 * Return the values of the load vector block
 * @param values: pointer to vector values
 * @return: false if load does not contribute
 *        vector element
 */
void Constantimpedance::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // load array starts from this location

  if(p_mode == RESIDUAL_EVAL) {
    f[0] = (p_va - p_R[0]*p_i[0])/p_L[0] - p_idot[0];
    f[1] = (p_vb - p_R[1]*p_i[1])/p_L[1] - p_idot[1];
    f[2] = (p_vc - p_R[2]*p_i[2])/p_L[2] - p_idot[2];
  }
}

  /**
   * Return the load current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
void Constantimpedance::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = p_i[0];
  *ib = p_i[1];
  *ic = p_i[2];
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Constantimpedance::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc;
}


/**
 * Get number of matrix values contributed by load
 * @return number of matrix values

 Non-zero pattern of the Jacobian is
         ia    ib    ic    va    vb    vc
 eq. 1 |  x                 x
 eq. 2 |        x                 x     
 eq. 3 |              x                 x

 Number of non-zeros in the Jacobian = 6
 */
int Constantimpedance::matrixNumValues()
{
  int numVals = 6;

  return numVals;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Constantimpedance::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  //partial w.r.t. load currents
  rows[ctr]   = p_gloc;
  rows[ctr+1] = p_gloc+1;
  rows[ctr+2] = p_gloc+2;
  
  cols[ctr]   = rows[ctr];
  cols[ctr+1] = rows[ctr+1];
  cols[ctr+2] = rows[ctr+2];

  values[ctr]   = -p_R[0]/p_L[0] - shift;
  values[ctr+1] = -p_R[1]/p_L[1] - shift;
  values[ctr+2] = -p_R[2]/p_L[2] - shift;

  ctr += 3;
  
  // Partial w.r.t voltages
  rows[ctr]   = p_gloc;
  rows[ctr+1] = p_gloc+1;
  rows[ctr+2] = p_gloc+2;

  cols[ctr]   = p_glocvoltage;
  cols[ctr+1] = p_glocvoltage+1;
  cols[ctr+2] = p_glocvoltage+2;

  values[ctr]   = 1.0/p_L[0];
  values[ctr+1] = 1.0/p_L[1];
  values[ctr+2] = 1.0/p_L[2];

  ctr += 3;

  *nvals = ctr;
  
}
