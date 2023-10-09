#include <gencls.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Gencls::Gencls(void)
{
  p_delta = 0.0;
  p_dw    = 0.0;
  p_deltadot = 0.0;
  p_dwdot    = 0.0;
  p_Rs    = 0.01;
  p_L     = 0.01/OMEGA_S;
  p_H     = 0.0;
  p_D     = 0.0;
  p_Ep    = 0.0;
  p_Pm    = 0.0;
  p_Xdp   = 0.0;
  
  nxgen   = 5; // Number of variables for this model
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

  // load parameters for the model type
  data->getValue(BUS_NUMBER, &bid);
  data->getValue(GENERATOR_RESISTANCE,&p_Rs,idx);
  data->getValue(GENERATOR_TRANSIENT_REACTANCE,&p_Xdp,idx);
  data->getValue(GENERATOR_INERTIA_CONSTANT_H,&p_H,idx);
  data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&p_D,idx);

  p_L = p_Xdp/OMEGA_S;
}

/**
 * Initialize generator model before calculation
 * @param [output] xin - array where initialized generator variables should be set
 */
void Gencls::init(gridpack::ComplexType* xin)
{
  double IGD,IGQ; // Machine currents in cartesian coordinates
  double Pg, Qg;  // Generator real and reactive power
  double dw=0.0;  // Initial machine speed deviation
  gridpack::ComplexType *x = xin+offsetb; // generator array starts from this location

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

  p_i[0] = ia;
  p_i[1] = ib;
  p_i[2] = ic;

  double p_Pm = (p_va*ia + p_vb*ib + p_vc*ic)*(2.0/3.0);
	
  x[0] = delta;
  x[1] = dw;
  x[2] = ia;
  x[3] = ib;
  x[4] = ic;
 
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
void Gencls::setValues(gridpack::ComplexType *values)
{
  gridpack::ComplexType *x = values+offsetb; // generator array starts from this location

  if(p_mode == XVECTOBUS) {
    p_delta = real(x[0]);
    p_dw    = real(x[1]);
    p_i[0]  = real(x[2]);
    p_i[1]  = real(x[3]);
    p_i[2]  = real(x[4]);
  } else if(p_mode == XDOTVECTOBUS) {
    p_deltadot = real(x[0]);
    p_dwdot    = real(x[1]);
    p_idot[0]  = real(x[2]);
    p_idot[1]  = real(x[3]);
    p_idot[2]  = real(x[4]);
  } 
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Gencls::vectorGetValues(gridpack::ComplexType *values)
{
  gridpack::ComplexType *f = values+offsetb; // generator array starts from this location

  if(p_mode == RESIDUAL_EVAL) {
    double Pe;
    double vabc[3];
    double vdq0[3];
    double idq0[3];
    double e[3];

    e[0] = p_Ep*sin(OMEGA_S*p_time + p_delta);
    e[1] = p_Ep*sin(OMEGA_S*p_time + p_delta - TWOPI_OVER_THREE);
    e[2] = p_Ep*sin(OMEGA_S*p_time + p_delta + TWOPI_OVER_THREE);

    vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;

    abc2dq0(vabc,p_time,p_delta,vdq0);
    abc2dq0(p_i,p_time,p_delta,idq0);

    Pe = vdq0[0]*idq0[0] + vdq0[1]*idq0[1] + vdq0[2]*idq0[2];
    
    // Generator equations
    f[0] = p_dw/OMEGA_S - p_deltadot;
    f[1]    = (p_Pm - Pe - p_D*p_dw)/(2*p_H) - p_dwdot;
    if(abs(p_L) > 1e-6) {
      // f = di_dt - idot => L^-1*(e - R*i - v) - idot
      f[2] = (e[0] - p_Rs*p_i[0] - p_va)/p_L - p_idot[0];
      f[3] = (e[1] - p_Rs*p_i[1] - p_vb)/p_L - p_idot[1];
      f[4] = (e[2] - p_Rs*p_i[2] - p_vc)/p_L - p_idot[2];
    } else {
      f[2] = (e[0] - p_Rs*p_i[0] - p_va);
      f[3] = (e[1] - p_Rs*p_i[1] - p_vb);
      f[4] = (e[2] - p_Rs*p_i[2] - p_vc);
    }      
    
    //printf("classical: %f\t%f\n", real(values[delta_idx]),real(values[dw_idx]));
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
  *ia = p_i[0];
  *ib = p_i[1];
  *ic = p_i[2];
}

/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values
 */
int Gencls::matrixNumValues()
{
}


/**
 * Return values from a matrix block
 * @param matrix - the Jacobian matrix
 */
void Gencls::matrixGetValues(gridpack::math::Matrix &matrix)
{
}


/**
 * Set Jacobian values
 * @param values a 2-d array of Jacobian block for the bus
 */
bool Gencls::setJacobian(gridpack::ComplexType **values)
{
  int VD_idx = 0; /* Row/col number for bus voltage VD variable */
  int VQ_idx = 1; /* Row/col number for bus voltage VQ variable */
  int IGQ_idx = 0; /* Row/col location for IGQ equations */
  int IGD_idx = 1; /* Row/col location for IGD equations */
  // Note that the order for the generator currents is [IGQ;IGD]
  int delta_idx = offsetb;
  int dw_idx = offsetb+1;

  /***** IMPORTANT NOTE ************/
  /********!!!!!!!!!!!!!!!!*********/
  // The array matrixDiagValues uses has a column-major-order! Because of
  // this we need to use the transpose of the locations. For e.g. the
  // partial derivative of dw equation w.r.t. VD would, in the natural
  // row-order form use a[dw_idx][VD_idx]. But, due to the column-major-
  // ordering, we need to use the transposed location, i.e., a[VD_idx][dw_idx]
  // This is extremely confusing!!!!!!!!

  if(p_mode == FAULT_EVAL) {
    // Generator variables held constant
    // dF_dX
    // Set diagonal values to 1.0
    values[delta_idx][delta_idx] = 1.0;
    values[dw_idx][dw_idx] = 1.0;

    // dG_dV
    values[VD_idx][IGQ_idx] +=  1/p_Xdp;
    values[VQ_idx][IGD_idx] += -1/p_Xdp;
  } else {

    // Partials of generator equations w.r.t generator variables
    // dF_dX
    values[delta_idx][delta_idx] = -shift;
    values[dw_idx][delta_idx] = 1.0/OMEGA_S;
    values[delta_idx][dw_idx] = (-VD*p_Ep*cos(p_delta)/p_Xdp - VQ*p_Ep*sin(p_delta)/p_Xdp)/(2*p_H);
    values[dw_idx][dw_idx] = -shift - p_D/(2*p_H);

    // dF_dV
    // These are the partial derivatives of the generator equations w.r.t voltage variables VD and VQ
    values[VD_idx][dw_idx] = (-p_Ep*sin(p_delta)/p_Xdp)/(2*p_H);
    values[VQ_idx][dw_idx] = (p_Ep*cos(p_delta)/p_Xdp)/(2*p_H);

    // dG_dX
    // These are the partial derivatives of the generator current w.r.t. generator variables
    values[delta_idx][IGQ_idx] = p_Ep*sin(p_delta)/p_Xdp;
    values[delta_idx][IGD_idx] = p_Ep*cos(p_delta)/p_Xdp;

    // dG_dV
    values[VD_idx][IGQ_idx] +=  1/p_Xdp;
    values[VQ_idx][IGD_idx] += -1/p_Xdp;
  }
			      
  return true;
}



