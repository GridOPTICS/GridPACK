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

BaseEMTExcModel::BaseEMTExcModel(void)
{
  p_nrows = 0;
  p_ncols = 0;
}

BaseEMTExcModel::~BaseEMTExcModel(void)
{
}

void BaseEMTExcModel::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx)
{
}

bool BaseEMTExcModel::setJacobian(gridpack::RealType **values)
{
  return false;
}

void BaseEMTExcModel::init(gridpack::RealType *values)
{
}

bool BaseEMTExcModel::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

void BaseEMTExcModel::write(const char* signal, char* string)
{
}

/**
 * Set Event
 */
void BaseEMTExcModel::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{

}

int BaseEMTExcModel::matrixNumValues()
{
  return 0;
}

void BaseEMTExcModel::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  *nvals = 0;
}

void BaseEMTExcModel::vectorGetValues(gridpack::RealType *values)
{
}

void BaseEMTExcModel::setValues(gridpack::RealType *values)
{
}

double BaseEMTExcModel::getFieldVoltage()
{
  return 0.0;
}

double BaseEMTExcModel::getFieldVoltage(int *Efd_gloc)
{
  *Efd_gloc = -1;
  return 0.0;
}


bool BaseEMTExcModel::getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen)
{
  return false;
}
