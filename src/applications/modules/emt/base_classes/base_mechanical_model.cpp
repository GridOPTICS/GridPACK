/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_mechanical_model.cpp
 *  
 * @brief Base Renewable generator mechanical part 
 *
 *
 */

#include <base_mechanical_model.hpp>
#include <gridpack/include/gridpack.hpp>

BaseEMTRMechModel::BaseEMTRMechModel(void)
{
  p_nrows = 0;
  p_ncols = 0;
}

BaseEMTRMechModel::~BaseEMTRMechModel(void)
{
}

void BaseEMTRMechModel::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx)
{
  data->getValue(BUS_NUMBER, &busnum);
  data->getValue(GENERATOR_ID, &id, idx); // Generator id
  data->getValue(GENERATOR_STAT,&status,idx); // Generator status
  data->getValue(CASE_SBASE,&sbase); // System MVAbase, used in conversion from machine base to system base.
  data->getValue(GENERATOR_MBASE,&mbase,idx); // Machine base (in MVA)
}

void BaseEMTRMechModel::init(gridpack::RealType *values)
{
}

bool BaseEMTRMechModel::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

void BaseEMTRMechModel::write(const char* signal, char* string)
{
}

/**
 * Set Event
 */
void BaseEMTRMechModel::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{

}

int BaseEMTRMechModel::matrixNumValues()
{
  return 0;
}

void BaseEMTRMechModel::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  *nvals = 0;
}

void BaseEMTRMechModel::vectorGetValues(gridpack::RealType *values)
{
}

void BaseEMTRMechModel::setValues(gridpack::RealType *values)
{
}
