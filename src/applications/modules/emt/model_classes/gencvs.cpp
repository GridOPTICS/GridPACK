#include <gencvs.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Gencvs::Gencvs(void)
{
  delta = 0.0;
  p_Rs    = 0.0;
  p_L     = 1.0/OMEGA_S;
  p_Ep    = 0.0;
  
  nxgen   = 3; // Number of variables for this model when integration type is implicit or implicit-explicit
}

void Gencvs::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgen = 3;
  *nvar = nxgen;
}

Gencvs::~Gencvs(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseGeneratorModel
 */
void Gencvs::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGenModel::load(data,idx); // load parameters in base generator model
  gridpack::ComplexType Zsource;

  // load parameters for the model type
  data->getValue(GENERATOR_ZSOURCE,&Zsource,idx);
  p_Rs = real(Zsource);
  p_Xdp = imag(Zsource);
  p_L = p_Xdp/OMEGA_S;
}

/**
 * Initialize generator model before calculation
 * @param [output] xin - array where initialized generator variables should be set
 */
void Gencvs::init(gridpack::RealType* xin)
{
  double IGD,IGQ; // Machine currents in cartesian coordinates
  double Pg, Qg;  // Generator real and reactive power
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

  delta = arg(E);
  
  double Im = abs(I);
  double Ia = arg(I);

  p_Ep = abs(E);

  double ia,ib,ic;

  ia = Im*sin(OMEGA_S*p_time + Ia);
  ib = Im*sin(OMEGA_S*p_time + Ia - TWOPI_OVER_THREE);
  ic = Im*sin(OMEGA_S*p_time + Ia + TWOPI_OVER_THREE);

  p_iabc[0] = ia;
  p_iabc[1] = ib;
  p_iabc[2] = ic;

  x[0] = ia;
  x[1] = ib;
  x[2] = ic;
}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Gencvs::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Gencvs::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void Gencvs::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values+offsetb; // generator array starts from this location

  if(p_mode == XVECTOBUS) {
    p_iabc[0]  = x[0];
    p_iabc[1]  = x[1];
    p_iabc[2]  = x[2];
  } else if(p_mode == XDOTVECTOBUS) {
    p_idot[0]  = x[0];
    p_idot[1]  = x[1];
    p_idot[2]  = x[2];
  } 
}

/**
   Prestep function
*/
void Gencvs::preStep(double time ,double timestep)
{

}

/**
   Poststep function
*/
void Gencvs::postStep(double time)
{
}



/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Gencvs::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // generator array starts from this location

  if(p_mode == RESIDUAL_EVAL) {
    double e[3];

    e[0] = p_Ep*sin(OMEGA_S*p_time + delta);
    e[1] = p_Ep*sin(OMEGA_S*p_time + delta - TWOPI_OVER_THREE);
    e[2] = p_Ep*sin(OMEGA_S*p_time + delta + TWOPI_OVER_THREE);

    p_vabc[0] = p_va;
    p_vabc[1] = p_vb;
    p_vabc[2] = p_vc;

    // Generator equations
    if(fabs(p_L) >= 1e-6) {
      // f = di_dt - idot => L^-1*(e - R*i - v) - idot
      f[0] = (e[0] - p_Rs*p_iabc[0] - p_va)/p_L - p_idot[0];
      f[1] = (e[1] - p_Rs*p_iabc[1] - p_vb)/p_L - p_idot[1];
      f[2] = (e[2] - p_Rs*p_iabc[2] - p_vc)/p_L - p_idot[2];
    } else {
      f[0] = (e[0] - p_Rs*p_iabc[0] - p_va);
      f[1] = (e[1] - p_Rs*p_iabc[1] - p_vb);
      f[2] = (e[2] - p_Rs*p_iabc[2] - p_vc);
    }
  }
}

  /**
   * Return the generator current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
void Gencvs::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = p_iabc[0];
  *ib = p_iabc[1];
  *ic = p_iabc[2];
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Gencvs::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc;
}


/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values

 Non-zero pattern of the Jacobian (x denotes non-zero value)

          ia   ib   ic    va    vb    vc
 eq. 1 |   x               x
 eq. 2 |        x               x
 eq. 3 |             x                x

 Number of non-zero values = 6
 */
int Gencvs::matrixNumValues()
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
void Gencvs::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  int ia_idx = p_gloc;
  int ib_idx = p_gloc+1;
  int ic_idx = p_gloc+2;

  int va_idx = p_glocvoltage;
  int vb_idx = p_glocvoltage+1;
  int vc_idx = p_glocvoltage+2;

  // partial derivatives w.r.t. f[0]-f[2]
  if(fabs(p_L) >= 1e-6) {
    rows[ctr]   = ia_idx;
    cols[ctr]   = ia_idx;
    values[ctr] = -p_Rs/p_L - shift;

    rows[ctr+1]   = ia_idx;
    cols[ctr+1]   = va_idx;
    values[ctr+1] = -1.0/p_L;

    ctr += 2;

    rows[ctr]   = ib_idx;
    cols[ctr]   = ib_idx;
    values[ctr] = -p_Rs/p_L - shift;

    rows[ctr+1]   = ib_idx;
    cols[ctr+1]   = vb_idx;
    values[ctr+1] = -1.0/p_L;

    ctr += 2;

    rows[ctr]   = ic_idx;
    cols[ctr]   = ic_idx;
    values[ctr] = -p_Rs/p_L - shift;

    rows[ctr+1]   = ic_idx;
    cols[ctr+1]   = vc_idx;
    values[ctr+1] = -1.0/p_L;

    ctr += 2;
  } else {
    rows[ctr]   = ia_idx;
    cols[ctr]   = ia_idx;
    values[ctr] = -p_Rs;

    rows[ctr+1]   = ia_idx;
    cols[ctr+1]   = va_idx;
    values[ctr+1] = -1.0;

    ctr += 2;

    rows[ctr]   = ib_idx;
    cols[ctr]   = ib_idx;
    values[ctr] = -p_Rs;

    rows[ctr+1]   = ib_idx;
    cols[ctr+1]   = vb_idx;
    values[ctr+1] = -1.0;

    ctr += 2;

    rows[ctr]   = ic_idx;
    cols[ctr]   = ic_idx;
    values[ctr] = -p_Rs;

    rows[ctr+1]   = ic_idx;
    cols[ctr+1]   = vc_idx;
    values[ctr+1] = -1.0;

    ctr += 2;
  }

  
  *nvals = ctr;
}
