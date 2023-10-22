#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_load_model.hpp"
#include <gridpack/include/gridpack.hpp>

BaseEMTLoadModel::BaseEMTLoadModel(void)
{
  pl = 0.0;
  ql = 0.0;
  mbase = DEFAULT_MVABASE;
  sbase = DEFAULT_MVABASE;
  status = 0;
  shift  = 0.0;
  p_va     = 0.0;
  p_vb     = 0.0;
  p_vc     = 0.0;
  nxload  = 0;
}

BaseEMTLoadModel::~BaseEMTLoadModel(void)
{
}

/**
 * Set Jacobian values
 * @param values a 2-d array of Jacobian block for the bus
 */
bool BaseEMTLoadModel::setJacobian(gridpack::ComplexType **values)
{
  return false;
}

/**
 * Load parameters from DataCollection object into load model
 * @param data collection of load parameters from input files
 * @param index of load on bus
 * TODO: might want to move this functionality to BaseEMTLoadModel
 */
void BaseEMTLoadModel::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx)
{
  data->getValue(BUS_NUMBER, &busnum);
  data->getValue(LOAD_STATUS,&status,idx); // Load status
  data->getValue(CASE_SBASE,&sbase); // System MVAbase, used in conversion from machine base to system base.
  if(status) {
    data->getValue(LOAD_PL,&pl,idx); // Load real power
    data->getValue(LOAD_QL,&ql,idx); // Load reactive power

    mbase = sbase; // Machine base is same as sbase    
  } else {
    pl = ql = mbase = 0.0;
  }
}

/**
 * Initialize load model before calculation
 * @param [output] values - array where initialized load variables should be set
 */
void BaseEMTLoadModel::init(gridpack::ComplexType *values)
{
}

/**
 * Write output from loads to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool BaseEMTLoadModel::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

/**
 * Write out load state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void BaseEMTLoadModel::write(const char* signal, char* string)
{
}

/**
 * Return the load current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void BaseEMTLoadModel::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = *ib = *ic = 0.0;
}

/**
 * Return the global location for the load current
 * @param [output] i_gloc - global location for the first current variable
 */
void BaseEMTLoadModel::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc;
}

/**
 * Get number of matrix values contributed by load
 * @return number of matrix values
 */
int BaseEMTLoadModel::matrixNumValues()
{
  return 0;
}

/**
 * Return values from a matrix block
 * @param matrix - the Jacobian matrix
 */
void BaseEMTLoadModel::matrixGetValues(gridpack::math::Matrix &matrix)
{
}

/**
 * Return vector values from the load model 
 * @param values - array of returned values
 *
 * Note: This function is used to return the entries in vector,
 * for e.g., the entries in the residual vector from the load
 * object
   */
void BaseEMTLoadModel::vectorGetValues(gridpack::ComplexType *values)
{
}

/**
 * Pass solution vector values to the load object
 * @param values - array of returned values
 *
 * Note: This function is used to pass the entries in vector
 * to the load object,
 * for e.g., the state vector values for this load
 */
void BaseEMTLoadModel::setValues(gridpack::ComplexType *values)
{
}


