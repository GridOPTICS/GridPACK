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

}

/**
 * Initialize generator model before calculation
 * @param [output] values - array where initialized generator variables should be set
 */
void Gencls::init(gridpack::ComplexType* values)
{
  double IGD,IGQ; // Machine currents in cartesian coordinates
  double Pg, Qg;  // Generator real and reactive power
  double delta,dw=0.0;  // Initial machine speed deviation
  double Vm;

  Pg = pg/sbase;
  Qg = qg/sbase;

  VD = Vm0*cos(Va0);
  VQ = Vm0*sin(Va0);

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
  if(p_mode == XVECTOBUS) {
    p_delta = real(values[0]);
    p_dw    = real(values[1]);
  } else if(p_mode == XDOTVECTOBUS) {
    p_deltadot = real(values[0]);
    p_dwdot    = real(values[1]);
  } else if(p_mode == XVECPRETOBUS) {
    p_deltaprev = real(values[0]);
    p_dwprev = real(values[1]);
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
  int delta_idx = 0, dw_idx = 1;
  // On fault (p_mode == FAULT_EVAL flag), the generator variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    //values[delta_idx] = values[dw_idx] = 0.0;
    values[delta_idx] = p_delta - p_deltaprev; 
    values[dw_idx] = p_dw - p_dwprev;
  } else if(p_mode == RESIDUAL_EVAL) {
    // Generator equations
    values[delta_idx] = p_dw/OMEGA_S - p_deltadot;
    values[dw_idx]    = (p_Pm - VD*p_Ep*sin(p_delta)/p_Xdp + VQ*p_Ep*cos(p_delta)/p_Xdp - p_D*p_dw)/(2*p_H) - p_dwdot;
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



