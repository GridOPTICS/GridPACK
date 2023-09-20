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
  p_nrows = 0;
  p_ncols = 0;
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
 * Set Jacobian block
 * @param value_map standard map containing indices and values of matrix
 *        elements
 */
bool BaseGenModel::setJacobian(std::map<std::pair<int,int>,gridpack::ComplexType>
    &value_map)
{
  return false;
}

#if 0
/**
 * Set the number of rows contributed by this generator
 * @param nrows number of rows
 */
void BaseGenModel::matrixSetNumRows(int nrows)
{
  p_nrows = nrows;
}

/**
 * Set the number of columns contributed by this generator
 * @param ncols number of columns
 */
void BaseGenModel::matrixSetNumCols(int ncols)
{
  p_ncols = ncols;
}

/**
 * Number of rows (equations) contributed to by generator
 * @return number of rows
 */
int BaseGenModel::matrixNumRows()
{
  return p_nrows;
}

/**
 * Number of rows (equations) contributed to by generator
 * @return number of rows
 */
int BaseGenModel::matrixNumCols()
{
  return p_ncols;
}

/** 
 * Set global row index
 * @param irow local row index
 * @param global row index
 */
void BaseGenModel::matrixSetRowIndex(int irow, int idx)
{
  if (p_rowidx.size() == 0) {
    p_rowidx.resize(p_nrows);
    int i;
    for (i=0; i<p_nrows; i++) p_rowidx[i] = -1;
  }
  p_rowidx[irow] = idx;
}

/** 
 * Set global column index
 * @param icol local column index
 * @param global column index
 */
void BaseGenModel::matrixSetColIndex(int icol, int idx)
{
  if (p_colidx.size() == 0) {
    p_colidx.resize(p_ncols);
    int i;
    for (i=0; i<p_ncols; i++) p_colidx[i] = -1;
  }
  p_colidx[icol] = idx;
}

/**
 * Return global row index given local row index
 * @param irow local row index
 * @return global row index
 */
int BaseGenModel::matrixGetRowIndex(int irow)
{
  return p_rowidx[irow];
}

/**
 * Return global column index given local column index
 * @param icol local column index
 * @return global column index
 */
int BaseGenModel::matrixGetColIndex(int icol)
{
  return p_colidx[icol];
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
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void BaseGenModel::matrixGetValues(int *nvals, gridpack::ComplexType *values,
    int *rows, int *cols)
{
}
#endif

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

