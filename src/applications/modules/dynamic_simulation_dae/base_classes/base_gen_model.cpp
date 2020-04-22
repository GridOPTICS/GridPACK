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
  VD     = 0.0;
  VQ     = 0.0;
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
  data->getValue(GENERATOR_STAT,&status,idx); // Generator status
  data->getValue(CASE_SBASE,&sbase); // System MVAbase, used in conversion from machine base to system base.
  if(status) {
    data->getValue(GENERATOR_PG,&pg,idx); // Generator real power
    data->getValue(GENERATOR_QG,&qg,idx); // Generator reactive power
    data->getValue(GENERATOR_MBASE,&mbase,idx); // Machine base (in MVA)
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

double BaseGenModel::getAngle()
{
  return 0.0;
}

/**
 *  Set the number of variables for this generator model
 *  @param [output] number of variables for this model
 */
bool BaseGenModel::vectorSize(int *nvar) const
{
  *nvar = 0;
  return true;
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
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void BaseGenModel::setValues(gridpack::ComplexType *values)
{
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
bool BaseGenModel::vectorValues(gridpack::ComplexType *values)
{
  return false;
}

/**
 * Return the generator current injection (in rectangular form) 
 * @param [output] IGD - real part of the generator current
 * @param [output] IGQ - imaginary part of the generator current
*/
void BaseGenModel::getCurrent(double *IGD, double *IGQ)
{
  *IGD = *IGQ = 0.0;
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

