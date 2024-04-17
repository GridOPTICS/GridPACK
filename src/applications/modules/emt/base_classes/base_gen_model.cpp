#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_gen_model.hpp"
#include <gridpack/include/gridpack.hpp>

BaseEMTGenModel::BaseEMTGenModel(void)
{
  pg = 0.0;
  qg = 0.0;
  mbase = DEFAULT_MVABASE;
  sbase = DEFAULT_MVABASE;
  status = 0;
  shift  = 0.0;
  p_va     = 0.0;
  p_vb     = 0.0;
  p_vc     = 0.0;
  nxgen  = 0;
  p_hasExciter = false;
  p_hasGovernor = false;
  p_hasPlantController = false;
}

BaseEMTGenModel::~BaseEMTGenModel(void)
{
}

/**
 * Set Jacobian values
 * @param values a 2-d array of Jacobian block for the bus
 */
bool BaseEMTGenModel::setJacobian(gridpack::RealType **values)
{
  return false;
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseEMTGeneratorModel
 */
void BaseEMTGenModel::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx)
{
  data->getValue(BUS_NUMBER, &busnum);
  data->getValue(GENERATOR_STAT,&status,idx); // Generator status
  data->getValue(CASE_SBASE,&sbase); // System MVAbase, used in conversion from machine base to system base.
  data->getValue(GENERATOR_MBASE,&mbase,idx); // Machine base (in MVA)
  if(status) {
    data->getValue(GENERATOR_PG,&pg,idx); // Generator real power
    data->getValue(GENERATOR_QG,&qg,idx); // Generator reactive power
    pg *= sbase;
    qg *= sbase;
  } else {
    pg = qg;
  }
}

/**
 * Initialize generator model before calculation
 * @param [output] values - array where initialized generator variables should be set
 */
void BaseEMTGenModel::init(gridpack::RealType *values)
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
bool BaseEMTGenModel::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void BaseEMTGenModel::write(const char* signal, char* string)
{
}

/**
 * Return the generator current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void BaseEMTGenModel::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = *ib = *ic = 0.0;
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void BaseEMTGenModel::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc;
}


/**
 * Get the field current
 * @param 
 */
double BaseEMTGenModel::getFieldCurrent()
{
  return 0.0;
}

/**
 * Return the rotor speed deviation
 * @param 
 */
double BaseEMTGenModel::getSpeedDeviation()
{
  return 0.0;
}

/**
 * Return the speed deviation and its global location 
 * @param[output] rotor speed deviation
 * @param[output] rotor speed deviation global location
 */
double BaseEMTGenModel::getSpeedDeviation(int *dw_gloc)
{
  *dw_gloc = 0;
  return 0.0;
}

void BaseEMTGenModel::setExciter(boost::shared_ptr<BaseEMTExcModel> &exciter)
{ 
  p_exciter = exciter;
  p_hasExciter = true;
}


boost::shared_ptr<BaseEMTExcModel> BaseEMTGenModel::getExciter()
{
  return p_exciter;
}

bool BaseEMTGenModel::hasExciter()
{
    return p_hasExciter;
}

void BaseEMTGenModel::setGovernor(boost::shared_ptr<BaseEMTGovModel> &governor)
{ 
  p_governor = governor;
  p_hasGovernor = true;
}

boost::shared_ptr<BaseEMTGovModel> BaseEMTGenModel::getGovernor()
{
  return p_governor;
}

bool BaseEMTGenModel::hasGovernor()
{
    return p_hasGovernor;
}

void BaseEMTGenModel::setPlantController(boost::shared_ptr<BaseEMTPlantControllerModel> &pcontroller)
{ 
  p_plantcontroller = pcontroller;
  p_hasPlantController = true;
}


boost::shared_ptr<BaseEMTPlantControllerModel> BaseEMTGenModel::getPlantController()
{
  return p_plantcontroller;
}

bool BaseEMTGenModel::hasPlantController()
{
    return p_hasPlantController;
}

/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values
 */
int BaseEMTGenModel::matrixNumValues()
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
void BaseEMTGenModel::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  *nvals = 0;
}

/**
 * Return vector values from the generator model 
 * @param values - array of returned values
 *
 * Note: This function is used to return the entries in vector,
 * for e.g., the entries in the residual vector from the generator
 * object
   */
void BaseEMTGenModel::vectorGetValues(gridpack::RealType *values)
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
void BaseEMTGenModel::setValues(gridpack::RealType *values)
{
}

/**
 * Returns the initial field voltage (Efd(t0))
 * @param [out] Efd0 - Initial field voltage
 */
double BaseEMTGenModel::getInitialFieldVoltage()
{
  return 0.0;
}

/**
   Prestep function
*/
void BaseEMTGenModel::preStep(double time, double timestep)
{
}

/**
   Poststep function
*/
void BaseEMTGenModel::postStep(double time)
{
}
