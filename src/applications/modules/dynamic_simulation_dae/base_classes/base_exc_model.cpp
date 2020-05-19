/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_exc_model.cpp
 * @author Shuangshuang Jin 
 * @author Shrirang Abhyankar
 * @Last modified:   04/22/20 - Shri
 *  
 * @brief Base exciter model 
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
}

/**
 * Set Jacobian block
 * @param values a 2-d array of Jacobian block for the bus
 */
bool BaseExcModel::setJacobian(gridpack::ComplexType **values)
{
  return false;
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

void BaseExcModel::setGenerator(BaseGenModel *generator)
{
  p_gen = generator;
}

BaseGenModel* BaseExcModel::getGenerator()
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

 
