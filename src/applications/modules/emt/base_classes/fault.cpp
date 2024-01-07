#include <fault.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Fault::Fault(void)
{
  nxfault   = 3; // Number of variables for this model

  Ron = 1e-2;
  Rgnd = 0.1;

  faulton = false;
}

Fault::~Fault(void)
{
}

void Fault::init(gridpack::ComplexType* xin)
{
  gridpack::ComplexType *x = xin + offsetb;

  x[0] = ifault[0] = 0.0;
  x[1] = ifault[1] = 0.0;
  x[2] = ifault[2] = 0.0;
}

/**
 * Write output from faults to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Fault::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out fault state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Fault::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto faults
 * @param values array containing fault state variables
*/
void Fault::setValues(gridpack::ComplexType *values)
{
  gridpack::ComplexType *x = values+offsetb; // load array starts from this location

  if(p_mode == XVECTOBUS) {
    ifault[0]  = real(x[0]);
    ifault[1]  = real(x[1]);
    ifault[2]  = real(x[2]);
  } 
}

/**
 * Return the values of the fault vector block
 * @param values: pointer to vector values
 * @return: false if fault does not contribute
 *        vector element
 */
void Fault::vectorGetValues(gridpack::ComplexType *values)
{
  gridpack::ComplexType *f = values+offsetb; // fault array starts from this location

  if(p_mode == RESIDUAL_EVAL) {
  }
}

  /**
   * Return the fault current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
void Fault::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = ifault[0];
  *ib = ifault[1];
  *ic = ifault[2];
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Fault::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc;
}


/**
 * Get number of matrix values contributed by fault
 * @return number of matrix values

 Non-zero pattern of the Jacobian is
         ia    ib    ic    va    vb    vc
 eq. 1 |  x                 x
 eq. 2 |        x                 x     
 eq. 3 |              x                 x

 Number of non-zeros in the Jacobian = 6
 */
int Fault::matrixNumValues()
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
void Fault::matrixGetValues(int *nvals, gridpack::ComplexType *values, int *rows, int *cols)
{
  int ctr = 0;

  //partial w.r.t. fault currents
  rows[ctr]   = p_gloc;
  rows[ctr+1] = p_gloc+1;
  rows[ctr+2] = p_gloc+2;
  
  cols[ctr]   = rows[ctr];
  cols[ctr+1] = rows[ctr+1];
  cols[ctr+2] = rows[ctr+2];

  //  values[ctr]   = -p_R[0];
  //  values[ctr+1] = -p_R[1];
  //  values[ctr+2] = -p_R[2];

  ctr += 3;
  
  // Partial w.r.t voltages
  rows[ctr]   = p_gloc;
  rows[ctr+1] = p_gloc+1;
  rows[ctr+2] = p_gloc+2;

  //  cols[ctr]   = p_glocvoltage;
  //  cols[ctr+1] = p_glocvoltage+1;
  //  cols[ctr+2] = p_glocvoltage+2;

  values[ctr]   = 1.0;
  values[ctr+1] = 1.0;
  values[ctr+2] = 1.0;

  ctr += 3;

  *nvals = ctr;
  
}


/**
 * Set Jacobian values
 * @param values a 2-d array of Jacobian block for the bus
 */
bool Fault::setJacobian(gridpack::ComplexType **values)
{
			      
  return true;
}

void Fault::setparams(EmtBus *emtbus,double faultontime, double faultofftime, std::string type, std::string phases, double Rfault, double Rg)
{
  bus = emtbus;
  ton = faultontime;
  toff = faultofftime;

  if(type == "SLG") {
    /* Single line to ground fault */
    faulttype = 1;
    faultedphases[0] = faultedphases[1] = faultedphases[2] = 0;
    if(phases == "A") faultedphases[0] = 1;
    if(phases == "B") faultedphases[1] = 1;
    if(phases == "C") faultedphases[2] = 1;

    Ron = Rfault;
    Rgnd = Rg;
  }
}




