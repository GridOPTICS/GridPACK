/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_gov_model.cpp
 *  
 * @brief  
 *
 *
 */

#include <base_gov_model.hpp>
#include <gridpack/include/gridpack.hpp>

BaseGovModel::BaseGovModel(void)
{
  mbase = DEFAULT_MVABASE;
  sbase = DEFAULT_MVABASE;
  p_nrows = 0;
  p_ncols = 0;
}

BaseGovModel::~BaseGovModel(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of goviteer on bus
 * TODO: might want to move this functionality to BaseGovModel
 */
void BaseGovModel::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx)
{
  data->getValue(CASE_SBASE,&sbase); // System MVAbase, used in conversion from machine base to system base.
  data->getValue(GENERATOR_MBASE,&mbase,idx); // Machine base (in MVA)
}

/**
 * Initialize governor model before calculation
 * @param [output] values - array where initialized governor variables should be set
 */
void BaseGovModel::init(gridpack::ComplexType *values)
{
}

/**
 * Write output from governors to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool BaseGovModel::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

/**
 *  Set the number of variables for this governor model
 *  @param [output] number of variables for this model
 */
bool BaseGovModel::vectorSize(int *nvar) const
{
  *nvar = 0;
  return true;
}

/**
 * Write out governor state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void BaseGovModel::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto governors
 * @param values array containing governor state variables
*/
void BaseGovModel::setValues(gridpack::ComplexType *values)
{
}

/**
 * Return the values of the governor vector block
 * @param values: pointer to vector values
 * @return: false if governor does not contribute
 *        vector element
 */
bool BaseGovModel::vectorValues(gridpack::ComplexType *values)
{
  return false;
}

/**
 * Set Jacobian block
 * @param values a 2-d array of Jacobian block for the bus
 */
bool BaseGovModel::setJacobian(gridpack::ComplexType **values)
{
  return false;
}

/**
 * Set Jacobian block
 * @param value_map standard map containing indices and values of matrix
 *        elements
 */
bool BaseGovModel::setJacobian(std::map<std::pair<int,int>,
    gridpack::ComplexType> &value_map)
{
  return false;
}

#if 0
/**
 * Set the number of rows contributed by this governor
 * @param nrows number of rows
 */
void BaseGovModel::matrixSetNumRows(int nrows)
{
  p_nrows = nrows;
}

/**
 * Set the number of columns contributed by this governor
 * @param ncols number of columns
 */
void BaseGovModel::matrixSetNumCols(int ncols)
{
  p_ncols = ncols;
}

/**
 * Number of rows (equations) contributed to by governor
 * @return number of rows
 */
int BaseGovModel::matrixNumRows()
{
  return p_nrows;
}

/**
 * Number of rows (equations) contributed to by governor
 * @return number of rows
 */
int BaseGovModel::matrixNumCols()
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
 *  * Return global row index given local row index
 *   * @param irow local row index
 *    * @return global row index
 *     */
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
 * Get number of matrix values contributed by governor
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
void BaseGenModel::matrixGetValues(int *nvals,gridpack::ComplexType *values,
    int *rows, int *cols)
{
}
#endif

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power
 */
void BaseGovModel::setInitialMechanicalPower(double pmech)
{
}

/** 
 * Get the value of the mechanical power parameter
 * @return value of the mechanical power
 */
double BaseGovModel::getMechanicalPower()
{
  return 0.0;
}

/**
 * Partial derivatives of Mechanical Power Pmech w.r.t. governor variables
 * @param xgov_loc locations of governor variables
 * @param dPmech_dxgov partial derivatives of mechanical power Pmech w.r.t governor variables
 */
bool BaseGovModel::getMechanicalPowerPartialDerivatives(int *xgov_loc,double *dPmech_dxgov)
{
  return false;
}


/** 
 * Set the value of VComp
 * @return value of Vcomp
 */
void BaseGovModel::setVcomp(double Vcomp)
{
}

/**
 * Set Event
 */
void BaseGovModel::setEvent(gridpack::math::DAESolver::EventManagerPtr eman)
{
}


void BaseGovModel::setGenerator(BaseGenModel *generator)
{
  p_gen = generator;
}

BaseGenModel* BaseGovModel::getGenerator(void)
{
  return p_gen;
}
 
