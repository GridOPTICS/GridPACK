#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_branch_model.hpp"
#include <gridpack/include/gridpack.hpp>

BaseEMTBranchModel::BaseEMTBranchModel(void)
{
  sbase = DEFAULT_MVABASE;
  status = 0;
  shift  = 0.0;
  nxbranch  = 0;
}

BaseEMTBranchModel::~BaseEMTBranchModel(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseEMTGeneratorModel
 */
void BaseEMTBranchModel::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx)
{
  // Get line
  data->getValue(BRANCH_STATUS,&status,idx);
  data->getValue(BRANCH_CKT,&cktid,idx);
  data->getValue(CASE_SBASE,&sbase); // System MVAbase
}

/**
 * Initialize branch model before calculation
 * @param [output] values - array where initialized branch variables should be set
 */
void BaseEMTBranchModel::init(gridpack::RealType *values)
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
bool BaseEMTBranchModel::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

/**
 * Write out branch state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void BaseEMTBranchModel::write(const char* signal, char* string)
{
}

/**
 * Return the branch from bus current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void BaseEMTBranchModel::getFromBusCurrent(double *ia, double *ib, double *ic)
{
  *ia = *ib = *ic = 0.0;
}

/**
 * Return the branch to bus current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void BaseEMTBranchModel::getToBusCurrent(double *ia, double *ib, double *ic)
{
  *ia = *ib = *ic = 0.0;
}

/**
 * Return the global location for the branch current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void BaseEMTBranchModel::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc;
}

/**
 * Get number of matrix values contributed by branch
 * @return number of matrix values
 */
int BaseEMTBranchModel::matrixNumValues()
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
void BaseEMTBranchModel::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
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
void BaseEMTBranchModel::vectorGetValues(gridpack::RealType *values)
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
void BaseEMTBranchModel::setValues(gridpack::RealType *values)
{
}

/**
   Prestep function
*/
void BaseEMTBranchModel::preStep(double time, double timestep)
{
}

/**
   Poststep function
*/
void BaseEMTBranchModel::postStep(double time)
{
}
