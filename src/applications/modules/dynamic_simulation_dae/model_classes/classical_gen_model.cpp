#include <classical_gen_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

ClassicalGen::ClassicalGen(void)
{
  p_delta = 0.0;
  p_dw    = 0.0;
  p_deltadot = 0.0;
  p_dwdot    = 0.0;
  p_Rs    = 0.0;
  p_H     = 0.0;
  p_D     = 0.0;
  p_Ep    = 0.0;
  p_Pm    = 0.0;
  p_Xdp   = 0.0;
}

ClassicalGen::~ClassicalGen(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseGeneratorModel
 */
void ClassicalGen::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseGenModel::load(data,idx); // load parameters in base generator model

  // load parameters for the model type
  data->getValue(GENERATOR_RESISTANCE,&p_Rs,idx);
  data->getValue(GENERATOR_TRANSIENT_REACTANCE,&p_Xdp,idx);
  data->getValue(GENERATOR_INERTIA_CONSTANT_H,&p_H,idx);
  data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&p_D,idx);

  // Convert generator parameters from machine base to MVA base
  p_H *= mbase/sbase;
  p_D *= mbase/sbase;
  p_Xdp /= mbase/sbase;

  printf("mbase = %f, sbase = %f\n", mbase, sbase);
}

/**
 * Initialize generator model before calculation
 * @param [output] values - array where initialized generator variables should be set
 */
void ClassicalGen::init(gridpack::ComplexType* values)
{
  double IGD,IGQ; // Machine currents in cartesian coordinates
  double Pg, Qg;  // Generator real and reactive power
  double delta,dw=0.0;  // Initial machine speed deviation
  double Vm;

  Pg = pg/sbase;
  Qg = qg/sbase;

  Vm = sqrt(VD*VD + VQ*VQ);

  IGD = (VD*Pg + VQ*Qg)/(Vm*Vm);
  IGQ = (VQ*Pg - VD*Qg)/(Vm*Vm);
  
  delta = atan2(VQ + p_Xdp*IGD,VD-p_Xdp*IGQ);

  p_Ep = sqrt(pow((VD - p_Xdp*IGQ),2) + pow((VQ + p_Xdp*IGD),2));
  p_Pm = Pg;
	
  values[0] = delta;
  values[1] = dw;

  printf("Pg = %f, Qg = %f, VD = %f, VQ = %f, Vm = %f, IGD = %f, IGQ = %f, delta = %f, p_Ep = %f, p_Pm = %f, dw = %f\n", Pg, Qg, VD, VQ, Vm, IGD, IGQ, delta, p_Ep, p_Pm, dw);
}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool ClassicalGen::serialWrite(char *string, const int bufsize,const char *signal)
{
}

double ClassicalGen::getAngle(void)
{
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void ClassicalGen::write(const char* signal, char* string)
{
}

/**
 *  Set the number of variables for this generator model
 *  @param [output] number of variables for this model
 */
bool ClassicalGen::vectorSize(int *nvar) const
{
  *nvar = 2;
  return true;
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void ClassicalGen::setValues(gridpack::ComplexType *values)
{
  if(p_mode == XVECTOBUS) {
    p_delta = real(values[0]);
    p_dw    = real(values[1]);
  } else if(p_mode == XDOTVECTOBUS) {
    p_deltadot = real(values[0]);
    p_dwdot    = real(values[1]);
  }
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
bool ClassicalGen::vectorValues(gridpack::ComplexType *values)
{
  int delta_idx = 0, dw_idx = 1;
  // On fault (p_mode == FAULT_EVAL flag), the generator variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    values[delta_idx] = values[dw_idx] = 0.0;
  } else if(p_mode == RESIDUAL_EVAL) {
    // Generator equations
    values[delta_idx] = p_dw/OMEGA_S - p_deltadot;
    values[dw_idx]    = (p_Pm - VD*p_Ep*sin(p_delta)/p_Xdp + VQ*p_Ep*cos(p_delta)/p_Xdp - p_D*p_dw)/(2*p_H) - p_dwdot;
    //printf("classical: %f\t%f\n", real(values[delta_idx]),real(values[dw_idx]));
  }
  
  return true;
}

/**
 * Return the generator current injection (in rectangular form) 
 * @param [output] IGD - real part of the generator current
 * @param [output] IGQ - imaginary part of the generator current
*/
void ClassicalGen::getCurrent(double *IGD, double *IGQ)
{
  // Generator current injections in the network
  *IGD += (-VQ + p_Ep*sin(p_delta))/p_Xdp;
  *IGQ += ( VD - p_Ep*cos(p_delta))/p_Xdp;
}

/**
 * Return the matrix entries
 * @param [output] nval - number of values set
 * @param [output] row - row indices for matrix entries
 * @param [output] col - col indices for matrix entries
 * @param [output] values - matrix entries
 * return true when matrix entries set
 */
bool ClassicalGen::matrixDiagEntries(int *nval,int *row, int *col, gridpack::ComplexType *values)
{
  int idx = 0;
  if(p_mode == FAULT_EVAL) {
    // On fault (p_mode == FAULT_EVAL flag), the generator variables are held constant. This is done by setting the diagonal matrix entries to 1.0 and all other entries to 0. The residual function values are already set to 0.0 in the vector values function. This results in the equation 1*dx = 0.0 such that dx = 0.0 and hence x does not get changed.
    row[idx] = 0; col[idx] = 0;
    values[idx] = 1.0;
    idx++;
    row[idx] = 1; col[idx] = 1;
    values[idx] = 1.0;
    idx++;
    *nval = idx;
  } else if(p_mode == DIG_DV) {
    // These are the partial derivatives of the generator currents (see getCurrent function) w.r.t to the voltage variables VD and VQ
    row[idx] = 0; col[idx] = 0;
    values[idx] =  1/p_Xdp;
    idx++;
    row[idx] = 1; col[idx] = 1;
    values[idx] =  -1/p_Xdp;
    idx++;

    *nval = idx;
  } else if(p_mode == DFG_DV) {
    // These are the partial derivatives of the generator equations w.r.t voltage variables VD and VQ
    row[idx] = 1; col[idx] = 0;
    values[idx] = (-p_Ep*sin(p_delta)/p_Xdp)/(2*p_H);
    idx++;
    row[idx] = 1; col[idx] = 1;
    values[idx] = (p_Ep*cos(p_delta)/p_Xdp)/(2*p_H);
    idx++;

    *nval = idx;
  } else if(p_mode == DIG_DX) {
    // These are the partial derivatives of the generator currents (see getCurrent) w.r.t. generator variables
    row[idx] = 0; col[idx] = 0;
    values[idx] = p_Ep*sin(p_delta)/p_Xdp;
    idx++;
    row[idx] = 1; col[idx] = 0;
    values[idx] = p_Ep*cos(p_delta)/p_Xdp;
    idx++;

    *nval = idx;
  } else {
    // Partials of generator equations w.r.t generator variables
    row[idx] = 0; col[idx] = 0;
    values[idx] = -shift;
    idx++;
    row[idx] = 0; col[idx] = 1;
    values[idx]    = 1.0/OMEGA_S;
    idx++;
    row[idx] = 1; col[idx] = 0;
    values[idx] = (-VD*p_Ep*cos(p_delta)/p_Xdp - VQ*p_Ep*sin(p_delta)/p_Xdp)/(2*p_H);
    idx++;
    row[idx] = col[idx] = 1;
    values[idx] = -shift - p_D/(2*p_H);
    idx++;
    
    *nval = idx;
  }
  return true;
}

