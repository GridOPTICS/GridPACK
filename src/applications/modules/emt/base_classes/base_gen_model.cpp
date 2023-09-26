#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_gen_model.hpp"
#include <gridpack/include/gridpack.hpp>

BaseGenModel::BaseGenModel(void)
{
  pg = 0.0;
  qg = 0.0;
  mbase = DEFAULT_MVABASE;
  sbase = DEFAULT_MVABASE;
  status = 0;
  shift  = 0.0;
  va     = 0.0;
  vb     = 0.0;
  vc     = 0.0;
  nxgen  = 0;
  p_hasExciter = false;
  p_hasGovernor = false;
}

BaseGenModel::~BaseGenModel(void)
{
}

/**
 * Set Jacobian values
 * @param values a 2-d array of Jacobian block for the bus
 */
bool BaseGenModel::setJacobian(gridpack::ComplexType **values)
{
  return false;
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseGeneratorModel
 */
void BaseGenModel::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx)
{
  data->getValue(BUS_NUMBER, &busnum);
  data->getValue(GENERATOR_STAT,&status,idx); // Generator status
  data->getValue(CASE_SBASE,&sbase); // System MVAbase, used in conversion from machine base to system base.
  if(status) {
    data->getValue(GENERATOR_PG,&pg,idx); // Generator real power
    data->getValue(GENERATOR_QG,&qg,idx); // Generator reactive power
    data->getValue(GENERATOR_MBASE,&mbase,idx); // Machine base (in MVA)
    pg *= sbase;
    qg *= sbase;
  } else {
    pg = qg = mbase = 0.0;
  }
}

/**
 * Initialize generator model before calculation
 * @param [output] values - array where initialized generator variables should be set
 */
void BaseGenModel::init(gridpack::ComplexType *values)
{
}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool BaseGenModel::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void BaseGenModel::write(const char* signal, char* string)
{
}

/**
 * Return the generator current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void BaseGenModel::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = *ib = *ic = 0.0;
}

/**
 * Get the field current
 * @param 
 */
double BaseGenModel::getFieldCurrent()
{
  return 0.0;
}

/**
 * Return the rotor speed deviation
 * @param 
 */
double BaseGenModel::getRotorSpeedDeviation()
{
  return 0.0;
}

/**
 * Return the location of speed rotor speed deviation variable in the bus array
 * @param rotor speed deviation location
*/
int BaseGenModel::getRotorSpeedDeviationLocation()
{
  return 0;
}

void BaseGenModel::setExciter(boost::shared_ptr<BaseExcModel> &exciter)
{ 
  p_exciter = exciter;
  p_hasExciter = true;
}


boost::shared_ptr<BaseExcModel> BaseGenModel::getExciter()
{
  return p_exciter;
}

bool BaseGenModel::hasExciter()
{
    return p_hasExciter;
}

void BaseGenModel::setGovernor(boost::shared_ptr<BaseGovModel> &governor)
{ 
  p_governor = governor;
  p_hasGovernor = true;
}

boost::shared_ptr<BaseGovModel> BaseGenModel::getGovernor()
{
  return p_governor;
}

bool BaseGenModel::hasGovernor()
{
    return p_hasGovernor;
}

/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values
 */
int BaseGenModel::matrixNumValues()
{
  return 0;
}

/**
 * Return values from a matrix block
 * @param matrix - the Jacobian matrix
 */
void BaseGenModel::matrixGetValues(gridpack::math::Matrix &matrix)
{
}

/**
 * Return vector values from the generator model 
 * @param values - array of returned values
 *
 * Note: This function is used to return the entries in vector,
 * for e.g., the entries in the residual vector from the generator
 * object
   */
void BaseGenModel::vectorGetValues(gridpack::ComplexType *values)
{
}

/**
 * Pass solution vector values to the generator object
 * @param values - array of returned values
 *
 * Note: This function is used to pass the entries in vector
 * to the generator object,
 * for e.g., the state vector values for this generator
 */
void BaseGenModel::setValues(gridpack::ComplexType *values)
{
}


