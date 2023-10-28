/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_exc_model.cpp
 *  
 * @brief Base exciter model 
 *
 *
 */

#include <base_exc_model.hpp>
#include <gridpack/include/gridpack.hpp>

BaseExcModel::BaseExcModel(void)
{
  p_nrows = 0;
  p_ncols = 0;
}

BaseExcModel::~BaseExcModel(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciteer on bus
 * TODO: might want to move this functionality to BaseExcModel
 */
void BaseExcModel::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx)
{
}


/**
 * Set the number of rows contributed by this excitor
 * @param nrows number of rows
 */
void BaseExcModel::matrixSetNumRows(int nrows)
{
  p_nrows = nrows;
}

/**
 * Set the number of columns contributed by this excitor
 * @param ncols number of columns
 */
void BaseExcModel::matrixSetNumCols(int ncols)
{
  p_ncols = ncols;
}

/**
 * Number of rows (equations) contributed to by excitor
 * @return number of rows
 */
int BaseExcModel::matrixNumRows()
{
  return p_nrows;
}

/**
 * Number of rows (equations) contributed to by excitor
 * @return number of rows
 */
int BaseExcModel::matrixNumCols()
{
  return p_ncols;
}

/** 
 * Set global row index
 * @param irow local row index
 * @param global row index
 */
void BaseExcModel::matrixSetRowIndex(int irow, int idx)
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
void BaseExcModel::matrixSetColIndex(int icol, int idx)
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
int BaseExcModel::matrixGetRowIndex(int irow)
{
    return p_rowidx[irow];
}

/**
 *  * Return global column index given local column index
 *   * @param icol local column index
 *    * @return global column index
 *     */
int BaseExcModel::matrixGetColIndex(int icol)
{
    return p_colidx[icol];
}

/**
 *  * Get number of matrix values contributed by excitor
 *   * @return number of matrix values
 *    */
int BaseExcModel::matrixNumValues()
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
void BaseExcModel::matrixGetValues(int *nvals,gridpack::RealType *values,
    int *rows, int *cols)
{
}

/**
 * Initialize exciter model before calculation
 * @param [output] values - array where initialized exciter variables should be set
 */
void BaseExcModel::init(gridpack::RealType *values)
{
}

/**
 * Write output from exciters to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool BaseExcModel::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

/**
 *  Set the number of variables for this exciter model
 *  @param [output] number of variables for this model
 */
bool BaseExcModel::vectorSize(int *nvar) const
{
  *nvar = 0;
  return true;
}

/**
 * Write out exciter state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void BaseExcModel::write(const char* signal, char* string)
{
}

/**
 * Set Event
 */
void BaseExcModel::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{

}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto exciters
 * @param values array containing exciter state variables
*/
void BaseExcModel::setValues(gridpack::RealType *values)
{
}

/**
 * Return the values of the exciter vector block
 * @param values: pointer to vector values
 * @return: false if exciter does not contribute
 *        vector element
 */
bool BaseExcModel::vectorValues(gridpack::RealType *values)
{
  return false;
}

/**
 * Set the initial field voltage (at t = tstart) parameter for the exciter
 * @param fldv value of the field voltage
 */
void BaseExcModel::setInitialFieldVoltage(double fldv)
{
}

/**
 * Partial derivatives of field voltage Efd w.r.t. exciter variables
 * @param xexc_loc locations of exciter variables
 * @param dEfd_dexc partial derivatives of field voltage w.r.t exciter variables
 * @param dEfd_dxgen partial derivatives of field voltage w.r.t. generator variables
 */
bool BaseExcModel::getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen)
{
  return false;
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double BaseExcModel::getFieldVoltage()
{
  return 0.0;
}

void BaseExcModel::setGenerator(BaseEMTGenModel *generator)
{
  p_gen = generator;
}

BaseEMTGenModel* BaseExcModel::getGenerator()
{
  return p_gen;
}


/**
 * Set the value of the time increment 
 * @return value of the time increment
 */
//void BaseExcModel::setTimeincrement(double timeincrement)
//{
//    t_inc = timeincrement;
//}

 
