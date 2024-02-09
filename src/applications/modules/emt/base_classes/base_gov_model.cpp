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

BaseEMTGovModel::BaseEMTGovModel(void)
{
  mbase = DEFAULT_MVABASE;
  sbase = DEFAULT_MVABASE;
  p_nrows = 0;
  p_ncols = 0;
}

BaseEMTGovModel::~BaseEMTGovModel(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of goviteer on bus
 * TODO: might want to move this functionality to BaseEMTGovModel
 */
void BaseEMTGovModel::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx)
{
  data->getValue(CASE_SBASE,&sbase); // System MVAbase, used in conversion from machine base to system base.
  data->getValue(GENERATOR_MBASE,&mbase,idx); // Machine base (in MVA)
}

/**
 * Initialize governor model before calculation
 * @param [output] values - array where initialized governor variables should be set
 */
void BaseEMTGovModel::init(gridpack::RealType *values)
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
bool BaseEMTGovModel::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

/**
 * Write out governor state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void BaseEMTGovModel::write(const char* signal, char* string)
{
}

/**
 * Set Jacobian block
 * @param values a 2-d array of Jacobian block for the bus
 */
bool BaseEMTGovModel::setJacobian(gridpack::RealType **values)
{
  return false;
}

/**
 * Return the values of the governor vector block
 * @param values: pointer to vector values
 * @return: false if governor does not contribute
 *        vector element
 */
void BaseEMTGovModel::vectorGetValues(gridpack::RealType *values)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto governors
 * @param values array containing governor state variables
*/
void BaseEMTGovModel::setValues(gridpack::RealType *values)
{
}

/**
 * Get number of matrix values contributed by governor
 * @return number of matrix values
 */
int BaseEMTGovModel::matrixNumValues()
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
void BaseEMTGovModel::matrixGetValues(int *nvals,gridpack::RealType *values,
    int *rows, int *cols)
{
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power
 */
void BaseEMTGovModel::setInitialMechanicalPower(double pmech)
{
}

/** 
 * Get the value of the mechanical power parameter
 * @return value of the mechanical power
 */
double BaseEMTGovModel::getMechanicalPower()
{
  return 0.0;
}

/** 
 * Get the value of the mechanical power and its global location
 * @return value of the mechanical power
 */
double BaseEMTGovModel::getMechanicalPower(int *Pmech_gloc)
{
  *Pmech_gloc = -1;
  return 0.0;
}


/**
 * Partial derivatives of Mechanical Power Pmech w.r.t. governor variables
 * @param xgov_loc locations of governor variables
 * @param dPmech_dxgov partial derivatives of mechanical power Pmech w.r.t governor variables
 */
bool BaseEMTGovModel::getMechanicalPowerPartialDerivatives(int *xgov_loc,double *dPmech_dxgov)
{
  return false;
}

/** 
 * Set the value of VComp
 * @return value of Vcomp
 */
void BaseEMTGovModel::setVcomp(double Vcomp)
{
}

/**
 * Set Event
 */
void BaseEMTGovModel::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
}
