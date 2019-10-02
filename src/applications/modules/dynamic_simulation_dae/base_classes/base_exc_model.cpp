/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_exc_model.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   08/23/19
 *  
 * @brief  
 *
 *
 */

#include <base_exc_model.hpp>
#include <gridpack/include/gridpack.hpp>

BaseExcModel::BaseExcModel(void)
{
  /*pg = 0.0;
  qg = 0.0;
  mbase = DEFAULT_MVABASE;
  sbase = DEFAULT_MVABASE;
  status = 0;
  shift  = 0.0;
  VD     = 0.0;
  VQ     = 0.0;
  hasExc = hasGovernor = false;*/
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
  /*data->getValue(GENERATOR_STAT,&status,idx); // Generator status
  data->getValue(CASE_SBASE,&sbase); // System MVAbase, used in conversion from machine base to system base.
  if(status) {
    data->getValue(GENERATOR_PG,&pg,idx); // Generator real power
    data->getValue(GENERATOR_QG,&qg,idx); // Generator reactive power
    data->getValue(GENERATOR_MBASE,&mbase,idx); // Machine base (in MVA)
  } else {
    pg = qg = mbase = 0.0;
  }*/
}

/**
 * Initialize exciter model before calculation
 * @param [output] values - array where initialized exciter variables should be set
 */
void BaseExcModel::init(gridpack::ComplexType *values)
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
}

double BaseExcModel::getAngle()
{
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
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto exciters
 * @param values array containing exciter state variables
*/
void BaseExcModel::setValues(gridpack::ComplexType *values)
{
}

/**
 * Return the values of the exciter vector block
 * @param values: pointer to vector values
 * @return: false if exciter does not contribute
 *        vector element
 */
bool BaseExcModel::vectorValues(gridpack::ComplexType *values)
{
  return false;
}

/**
 * Return the exciter current injection (in rectangular form) 
 * @param [output] IGD - real part of the exciter current
 * @param [output] IGQ - imaginary part of the exciter current
*/
void BaseExcModel::getCurrent(double *IGD, double *IGQ)
{
  *IGD = *IGQ = 0.0;
}

/**
 * Return the matrix entries
 * @param [output] nval - number of values set
 * @param [output] row - row indices for matrix entries
 * @param [output] col - col indices for matrix entries
 * @param [output] values - matrix entries
 * return true when matrix entries set
 */
bool BaseExcModel::matrixDiagEntries(int *nval,int *row, int *col, gridpack::ComplexType *values)
{
  *nval = 0;
  return false;
}

/**
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void BaseExcModel::setFieldVoltage(double fldv)
{
}

/**
 * Set the field current parameter inside the exciter
 * @param fldc value of the field current
 */
void BaseExcModel::setFieldCurrent(double fldc)
{
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double BaseExcModel::getFieldVoltage()
{
  return 0.0;
}

/** 
 * Get the value of the field current parameter
 * @return value of field current
 */
double BaseExcModel::getFieldCurrent()
{
  return 0.0;
}

/** 
 * Set the value of the Vterminal
 * @return value of field current
 */
void BaseExcModel::setVterminal(double mag)
{
}

/** 
 * Set the value of the omega
 * @return value of field current
 */
void BaseExcModel::setOmega(double omega)
{
}

/** 
 * Set the value of VComp
 * @return value of Vcomp
 */
void BaseExcModel::setVcomp(double Vcomp)
{
}

/**
 * Set the value of the time step
 * @return value of the time step
 */
void BaseExcModel::setTimestep(double timestep)
{
    ts = timestep;
}

/**
 * Set the value of the time increment 
 * @return value of the time increment
 */
void BaseExcModel::setTimeincrement(double timeincrement)
{
    t_inc = timeincrement;
}

 
