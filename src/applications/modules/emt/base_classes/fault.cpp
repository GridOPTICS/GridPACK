#include <fault.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Fault::Fault(void)
{
  nxfault   = 3; // Number of variables for this model

  Ron = 1e-2;
  Rgnd = 0.1;

  faulton[0] = faulton[1] = faulton[2] = false;
  faultoff = false;
}

Fault::~Fault(void)
{
}

void Fault::init(gridpack::RealType* xin)
{
  gridpack::RealType *x = xin + offsetb;

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
void Fault::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values+offsetb; // load array starts from this location

  if(p_mode == XVECTOBUS) {
    ifault[0]  = x[0];
    ifault[1]  = x[1];
    ifault[2]  = x[2];
  } 
}

/**
 * Return the values of the fault vector block
 * @param values: pointer to vector values
 * @return: false if fault does not contribute
 *        vector element
 */
void Fault::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // fault array starts from this location

  if(p_mode == RESIDUAL_EVAL) {
    f[0] = ifault[0];
    f[1] = ifault[1];
    f[2] = ifault[2];

    if(faulttype == 1) { // SLG
      if(faulton[0] && faultedphases[0]) f[0] = p_va - (Ron + Rgnd)*ifault[0];
      if(faulton[1] && faultedphases[1]) f[1] = p_vb - (Ron + Rgnd)*ifault[1];
      if(faulton[2] && faultedphases[2]) f[2] = p_vc - (Ron + Rgnd)*ifault[2];

    } else if(faulttype == 3) { // ThreePhase
      if(faulton[0]) f[0] = p_va - (Ron + Rgnd)*ifault[0] - Rgnd*ifault[1]         - Rgnd*ifault[2];
      if(faulton[1]) f[1] = p_vb -         Rgnd*ifault[0] - (Ron + Rgnd)*ifault[1] - Rgnd*ifault[2];
      if(faulton[2]) f[2] = p_vc -         Rgnd*ifault[0] - Rgnd*ifault[1]         - (Ron + Rgnd)*ifault[2];	
    }
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
  int numVals = 0;
  if(faulttype == 1) {
    numVals = 6;
  } else if (faulttype == 3) {
    numVals = 12;
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
void Fault::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  if(faulttype == 1) {
    //partial w.r.t. fault currents
    rows[ctr]   = p_gloc;
    rows[ctr+1] = p_gloc+1;
    rows[ctr+2] = p_gloc+2;
    
    cols[ctr]   = rows[ctr];
    cols[ctr+1] = rows[ctr+1];
    cols[ctr+2] = rows[ctr+2];

    values[ctr] = 1.0;
    values[ctr+1] = 1.0;
    values[ctr+2] = 1.0;
  } else if(faulttype == 3) {
    rows[ctr]   = p_gloc;      cols[ctr]   = p_gloc;
    rows[ctr+1] = p_gloc;      cols[ctr+1] = p_gloc + 1;
    rows[ctr+2] = p_gloc;      cols[ctr+2] = p_gloc + 2;
    rows[ctr+3] = p_gloc + 1;  cols[ctr+3] = p_gloc;
    rows[ctr+4] = p_gloc + 1;  cols[ctr+4] = p_gloc + 1;
    rows[ctr+5] = p_gloc + 1;  cols[ctr+5] = p_gloc + 2;
    rows[ctr+6] = p_gloc + 2;  cols[ctr+6] = p_gloc;
    rows[ctr+7] = p_gloc + 2;  cols[ctr+7] = p_gloc + 1;
    rows[ctr+8] = p_gloc + 2;  cols[ctr+8] = p_gloc + 2;

    values[ctr]   = 1.0;
    values[ctr+1] = 0.0;
    values[ctr+2] = 0.0;
    values[ctr+3] = 0.0;
    values[ctr+4] = 1.0;
    values[ctr+5] = 0.0;
    values[ctr+6] = 0.0;
    values[ctr+7] = 0.0;
    values[ctr+8] = 1.0;
  }

  if(faulttype == 1) { // SLG
    if(faulton[0] && faultedphases[0]) values[ctr] = -(Ron + Rgnd);
    if(faulton[1] && faultedphases[1]) values[ctr+1] = -(Ron + Rgnd);
    if(faulton[2] && faultedphases[2]) values[ctr+2] = -(Ron + Rgnd);

    ctr += 3;
  } else if(faulttype == 3) { //ThreePhase{
    if(faulton[0]) {
      values[ctr]   = -(Ron + Rgnd);
      values[ctr+1] = -Rgnd;
      values[ctr+2] = -Rgnd;
    }

    if(faulton[1]) {
      values[ctr+3] = -Rgnd;
      values[ctr+4] = -(Ron + Rgnd);;
      values[ctr+5] = -Rgnd;
    }

    if(faulton[2]) {
      values[ctr+6] = -Rgnd;
      values[ctr+7] = -Rgnd;
      values[ctr+8] = -(Ron + Rgnd);
    }
    ctr += 9;
  }

  rows[ctr]   = p_gloc;
  rows[ctr+1] = p_gloc+1;
  rows[ctr+2] = p_gloc+2;
  
  cols[ctr]   = p_glocvoltage;
  cols[ctr+1] = p_glocvoltage+1;
  cols[ctr+2] = p_glocvoltage+2;
  
  values[ctr] = values[ctr+1] = values[ctr+2] = 0.0;

  if(faulttype == 1) { // SLG
    // Partial w.r.t voltages
    if(faulton[0] && faultedphases[0]) values[ctr]   = 1.0;
    if(faulton[1] && faultedphases[1]) values[ctr+1] = 1.0;
    if(faulton[2] && faultedphases[2]) values[ctr+2] = 1.0;
  } else if(faulttype == 3) { // Three Phase
    if(faulton[0]) values[ctr] = 1.0;
    if(faulton[1]) values[ctr+1] = 1.0;
    if(faulton[2]) values[ctr+2] = 1.0;
  }

  ctr += 3;

  *nvals = ctr;
  
}


/**
 * Set Jacobian values
 * @param values a 2-d array of Jacobian block for the bus
 */
bool Fault::setJacobian(gridpack::RealType **values)
{
			      
  return true;
}

void Fault::setparams(double faultontime, double faultofftime, std::string type, std::string phases, double Rfault, double Rg)
{
  ton = faultontime;
  toff = faultofftime;

  if(type == "SLG") {
    /* Single line to ground fault */
    faulttype = 1;
    faultedphases[0] = faultedphases[1] = faultedphases[2] = 0;
    if(phases == "A") faultedphases[0] = 1;
    if(phases == "B") faultedphases[1] = 1;
    if(phases == "C") faultedphases[2] = 1;
  } else if(type == "ThreePhase") {
    /* Three phase fault */
    faulttype = 3;
    faultedphases[0] = faultedphases[1] = faultedphases[2] = 1;
  }
  
  Ron = Rfault;
  Rgnd = Rg;
  
}

void Fault::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  gridpack::math::RealDAESolver::EventPtr e(new FaultEvent(this));

  eman->add(e);

}

/**
 * Update the event function values
 */
void Fault::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
  int offset    = getLocalOffset();

  ifault[0] = state[offset];
  ifault[1] = state[offset+1];
  ifault[2] = state[offset+2];

  evalues[0] = ton - t;
  evalues[1] = toff - t;

  evalues[2] = evalues[3] = evalues[4] = 1.0;
  if(faultoff) { /* Time greater than fault off time */
    if(faulton[0] && faultedphases[0]) evalues[2] = ifault[0];
    if(faulton[1] && faultedphases[1]) evalues[3] = ifault[1];
    if(faulton[2] && faultedphases[2]) evalues[4] = ifault[2];
  }
} 

/**
 * Event handler
 */
void Fault::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int offset    = getLocalOffset();

  ifault[0] = state[offset];
  ifault[1] = state[offset+1];
  ifault[2] = state[offset+2];

  if(triggered[0]) { /* Fault on */
    if(faulttype == 1) {
      if(faultedphases[0]) faulton[0] = true;
      if(faultedphases[1]) faulton[1] = true;
      if(faultedphases[2]) faulton[2] = true;
    } else if(faulttype == 3) {
      faulton[0] = faulton[1] = faulton[2] = true;
    }
  }

  if(triggered[1]) { /* Fault is off */
    faultoff = true;
  }

  if(triggered[2]) {
    faulton[0] = false; /* Phase A fault current reached zero */
  }
  if(triggered[3]) {
    faulton[1] = false; /* Phase B fault current reached zero */
  }
  if(triggered[4]) {
    faulton[2] = false; /* Phase C fault current reached zero */
  }

  if(faultoff) {
    if((faulton[0] + faulton[1] + faulton[2]) == false) { /* Fault has extinguished */
      faultoff = false;
    }
  }

  
}

void FaultEvent::p_update(const double& t,gridpack::RealType *state)
{
  p_fault->eventFunction(t,state,p_current);
}

void FaultEvent::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_fault->eventHandlerFunction(triggered,t,state);
}

