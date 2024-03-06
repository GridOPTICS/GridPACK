#include <lumpedline.hpp>
#include <gridpack/include/gridpack.hpp>

Lumpedline::Lumpedline(void)
{
  nxbranch  = 3;
}

Lumpedline::~Lumpedline(void)
{
}

/**
 * Return the number of variables
 * @param [output] nvar - number of variables
 */
void Lumpedline::getnvar(int *nvar)
{
  *nvar = nxbranch;
}


/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseEMTGeneratorModel
 */
void Lumpedline::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int i)
{
  // Get line
  BaseEMTBranchModel::load(data,i);

  data->getValue(BRANCH_R,&R,i);
  data->getValue(BRANCH_X,&X,i);
  data->getValue(BRANCH_B,&Bc,i);

  if(abs(R) > 1e-6) {
    p_hasResistance = true;
  }

  if(abs(X) > 1e-6) {
    p_hasInductance = true;
  }

  R1 = R; L1 = X/OMEGA_S; C1 = Bc/OMEGA_S;
  R0 = 3*R1; L0 = 3*L1; C0 = 3*C1;

  double Rs = (2*R1 + R0)/3.0;
  double Rm = (R0 - R1)/3.0;
  p_R[0][0] = p_R[1][1] = p_R[2][2] = Rs;
  p_R[0][1] = p_R[1][0] = Rm;
  p_R[0][2] = p_R[2][0] = Rm;
  p_R[1][2] = p_R[2][1] = Rm;
  
  double Ls = (2*L1 + L0)/3.0;
  double Lm = (L0 - L1)/3.0;
  p_L[0][0] = p_L[1][1] = p_L[2][2] = Ls;
  p_L[0][1] = p_L[1][0] = Lm;
  p_L[0][2] = p_L[2][0] = Lm;
  p_L[1][2] = p_L[2][1] = Lm;
  
  double Cp = (2*C1 + C0)/3.0;
  double Cg = (C0 - C1)/3.0;
  p_C[0][0] = p_C[1][1] = p_C[2][2] = Cp;
  p_C[0][1] = p_C[1][0] = Cg;
  p_C[0][2] = p_C[2][0] = Cg;
  p_C[1][2] = p_C[2][1] = Cg;
}

/**
 * Initialize branch model before calculation
 * @param [output] values - array where initialized branch variables should be set
 */
void Lumpedline::init(gridpack::RealType *values)
{
}

/**
 * Write output from branchs to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Lumpedline::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

/**
 * Write out branch state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Lumpedline::write(const char* signal, char* string)
{
}

/**
 * Return the branch current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void Lumpedline::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = *ib = *ic = 0.0;
}

/**
 * Return the global location for the branch current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Lumpedline::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc;
}

/**
 * Get number of matrix values contributed by branch
 * @return number of matrix values
 */
int Lumpedline::matrixNumValues()
{
  return 0;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Lumpedline::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  *nvals = 0;
}

/**
 * Return vector values from the branch model 
 * @param values - array of returned values
 *
 * Note: This function is used to return the entries in vector,
 * for e.g., the entries in the residual vector from the branch
 * object
   */
void Lumpedline::vectorGetValues(gridpack::RealType *values)
{
}

/**
 * Pass solution vector values to the branch object
 * @param values - array of returned values
 *
 * Note: This function is used to pass the entries in vector
 * to the branch object,
 * for e.g., the state vector values for this branch
 */
void Lumpedline::setValues(gridpack::RealType *values)
{
}

/**
   Prestep function
*/
void Lumpedline::preStep(double time, double timestep)
{
}

/**
   Poststep function
*/
void Lumpedline::postStep(double time)
{
}
