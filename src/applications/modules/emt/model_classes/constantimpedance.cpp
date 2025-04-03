#include <constantimpedance.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Constantimpedance::Constantimpedance(void)
{
  series_RL = false;

  if(series_RL) {
    nxload   = 3; // Number of variables for this model (assuming the reactive part of load exists)
  } else {
    nxload = 6;
  }
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
  double Im,Ia,ILm,ILa;
  
  gridpack::RealType *x = xin + offsetb;

  VD = p_Vm0*cos(p_Va0);
  VQ = p_Vm0*sin(p_Va0);

  gridpack::ComplexType V = gridpack::ComplexType(VD,VQ);
  gridpack::ComplexType S = gridpack::ComplexType(pl,ql);
  gridpack::ComplexType I = conj(S/V);
  gridpack::ComplexType Z = V/I;
  gridpack::ComplexType Y = I/V;

  if(series_RL) {
    R = real(Z);
    L = imag(Z)/OMEGA_S;

    Im = abs(I);
    Ia = arg(I);
  } else {

    double ZR,ZI,ZI_ZR;
    double XL;

    ZR = real(Z);
    ZI = imag(Z);

    ZI_ZR = ZI/ZR;
    
    R = (ZR*ZR + ZI*ZI)/ZR;
    XL = (ZR*ZR + ZI*ZI)/ZI;

    L = XL/OMEGA_S;

    // Parallel RL
    //    R = 1/real(Y);
    //    L = -1/imag(Y)/OMEGA_S;

    gridpack::ComplexType IL;
    gridpack::ComplexType Jay = gridpack::ComplexType(0.0,1.0);
    IL = I - V/R;

    ILm = abs(IL);
    ILa = arg(IL);
  }    


  p_R[0] = p_R[1] = p_R[2] = R;
  p_L[0] = p_L[1] = p_L[2] = L;

  if(series_RL) {
    x[0] = p_i[0] = Im*sin(Ia);
    x[1] = p_i[1] = Im*sin(Ia - TWOPI_OVER_THREE);
    x[2] = p_i[2] = Im*sin(Ia + TWOPI_OVER_THREE);
  } else {
    x[0] = p_i[0] = ILm*sin(ILa);
    x[1] = p_i[1] = ILm*sin(ILa - TWOPI_OVER_THREE);
    x[2] = p_i[2] = ILm*sin(ILa + TWOPI_OVER_THREE);

    Im = abs(I);
    Ia = arg(I);
    
    x[3] = p_i[3] = Im*sin(Ia);
    x[4] = p_i[4] = Im*sin(Ia - TWOPI_OVER_THREE);
    x[5] = p_i[5] = Im*sin(Ia + TWOPI_OVER_THREE);
  }
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

  if(series_RL) {
    if(p_mode == XVECTOBUS) {
      p_i[0]  = x[0];
      p_i[1]  = x[1];
      p_i[2]  = x[2];
    } else if(p_mode == XDOTVECTOBUS) {
      p_idot[0]  = x[0];
      p_idot[1]  = x[1];
      p_idot[2]  = x[2];
    }
  } else {
    if(p_mode == XVECTOBUS) {
      p_i[0]  = x[0];
      p_i[1]  = x[1];
      p_i[2]  = x[2];
      p_i[3]  = x[3];
      p_i[4]  = x[4];
      p_i[5]  = x[5];
    } else if(p_mode == XDOTVECTOBUS) {
      p_idot[0]  = x[0];
      p_idot[1]  = x[1];
      p_idot[2]  = x[2];
    }
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
    if(series_RL) {
      if(fabs(ql) <= 1e-6) {// no ql
	f[0] = (p_va - p_R[0]*p_i[0]);
	f[1] = (p_vb - p_R[1]*p_i[1]);
	f[2] = (p_vc - p_R[2]*p_i[2]);
      } else {
	f[0] = (p_va - p_R[0]*p_i[0])/p_L[0] - p_idot[0];
	f[1] = (p_vb - p_R[1]*p_i[1])/p_L[1] - p_idot[1];
	f[2] = (p_vc - p_R[2]*p_i[2])/p_L[2] - p_idot[2];
      }
    } else {
      f[0] = p_va/p_L[0] - p_idot[0];
      f[1] = p_vb/p_L[1] - p_idot[1];
      f[2] = p_vc/p_L[2] - p_idot[2];

      f[3] = p_i[0] + p_va/p_R[0] - p_i[3];
      f[4] = p_i[1] + p_vb/p_R[1] - p_i[4];
      f[5] = p_i[2] + p_vc/p_R[2] - p_i[5];
    }
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
  if(series_RL) {
    *ia = p_i[0];
    *ib = p_i[1];
    *ic = p_i[2];
  } else {
    *ia = p_i[3];
    *ib = p_i[4];
    *ic = p_i[5];
  }
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Constantimpedance::getCurrentGlobalLocation(int *i_gloc)
{
  if(series_RL) {
    *i_gloc = p_gloc;
  } else {
    *i_gloc = p_gloc + 3;
  }
}


/**
 * Get number of matrix values contributed by load
 * @return number of matrix values

Series_RL
---------
 Non-zero pattern of the Jacobian is
         ia    ib    ic    va    vb    vc
 eq. 1 |  x                 x
 eq. 2 |        x                 x     
 eq. 3 |              x                 x

 Number of non-zeros in the Jacobian = 6

Parallel RL
-----------
  Non-zero pattern of the Jacobian is
         iLa    iLb    iLc   ia   ib   ic    va    vb    vc
 eq. 1 |  x                                   x
 eq. 2 |         x                                  x     
 eq. 3 |                x                                 x
 eq. 4 |  x                   x               x
 eq. 5 |         x                   x               x
 eq. 6 |                x                x                x
  Number of non-zeros in the Jacobian = 15
 

 */
int Constantimpedance::matrixNumValues()
{
  int numVals;
  if(series_RL) {
    numVals = 6;
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
void Constantimpedance::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  if(series_RL) {
    //partial w.r.t. load currents
    rows[ctr]   = p_gloc;
    rows[ctr+1] = p_gloc+1;
    rows[ctr+2] = p_gloc+2;
    
    cols[ctr]   = rows[ctr];
    cols[ctr+1] = rows[ctr+1];
    cols[ctr+2] = rows[ctr+2];
    
    if(fabs(ql) <= 1e-6) {
      values[ctr]   = -p_R[0];
      values[ctr+1] = -p_R[1];
      values[ctr+2] = -p_R[2];
    } else {
      values[ctr]   = -p_R[0]/p_L[0] - shift;
      values[ctr+1] = -p_R[1]/p_L[1] - shift;
      values[ctr+2] = -p_R[2]/p_L[2] - shift;
    }
    
    ctr += 3;
    
    // Partial w.r.t voltages
    rows[ctr]   = p_gloc;
    rows[ctr+1] = p_gloc+1;
    rows[ctr+2] = p_gloc+2;
    
    cols[ctr]   = p_glocvoltage;
    cols[ctr+1] = p_glocvoltage+1;
    cols[ctr+2] = p_glocvoltage+2;
    
    if(fabs(ql) <= 1e-6) {
      values[ctr]   = 1.0;
      values[ctr+1] = 1.0;
      values[ctr+2] = 1.0;
    } else {
      values[ctr]   = 1.0/p_L[0];
      values[ctr+1] = 1.0/p_L[1];
      values[ctr+2] = 1.0/p_L[2];
    }
    ctr += 3;
  } else { // Parallel RL
    //partial eq 0-2 w.r.t. inductor currents
    rows[ctr]   = p_gloc;
    rows[ctr+1] = p_gloc+1;
    rows[ctr+2] = p_gloc+2;
    
    cols[ctr]   = p_gloc;
    cols[ctr+1] = p_gloc+1;
    cols[ctr+2] = p_gloc+2;
    
    values[ctr]   = -shift;
    values[ctr+1] = -shift;
    values[ctr+2] = -shift;
    
    ctr += 3;

    // Partial eq 0-2 w.r.t voltages
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

    //partial eq 3-5 w.r.t. inductor currents
    rows[ctr]   = p_gloc+3;
    rows[ctr+1] = p_gloc+4;
    rows[ctr+2] = p_gloc+5;
    
    cols[ctr]   = p_gloc;
    cols[ctr+1] = p_gloc+1;
    cols[ctr+2] = p_gloc+2;
    
    values[ctr]   = 1.0;
    values[ctr+1] = 1.0;
    values[ctr+2] = 1.0;
    
    ctr += 3;

    // Partial eq 3-5 w.r.t voltages
    rows[ctr]   = p_gloc+3;
    rows[ctr+1] = p_gloc+4;
    rows[ctr+2] = p_gloc+5;
    
    cols[ctr]   = p_glocvoltage;
    cols[ctr+1] = p_glocvoltage+1;
    cols[ctr+2] = p_glocvoltage+2;
    
    values[ctr]   = 1.0/p_R[0];
    values[ctr+1] = 1.0/p_R[1];
    values[ctr+2] = 1.0/p_R[2];
    
    ctr += 3;

    //partial eq 3-5 w.r.t. load currents
    rows[ctr]   = p_gloc+3;
    rows[ctr+1] = p_gloc+4;
    rows[ctr+2] = p_gloc+5;
    
    cols[ctr]   = p_gloc+3;
    cols[ctr+1] = p_gloc+4;
    cols[ctr+2] = p_gloc+5;
    
    values[ctr]   = -1.0;
    values[ctr+1] = -1.0;
    values[ctr+2] = -1.0;
    
    ctr += 3;
    
  }
    
  *nvals = ctr;
  
}


/**
 * Set Jacobian values
 * @param values a 2-d array of Jacobian block for the bus
 */
bool Constantimpedance::setJacobian(gridpack::RealType **values)
{
			      
  return true;
}



