/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_plant_model.cpp
 *  
 * @brief Base plant model 
 *
 *
 */

#include <base_plant_model.hpp>
#include <gridpack/include/gridpack.hpp>

BaseEMTPlantControllerModel::BaseEMTPlantControllerModel(void)
{
  p_nrows = 0;
  p_ncols = 0;
}

BaseEMTPlantControllerModel::~BaseEMTPlantControllerModel(void)
{
}

void BaseEMTPlantControllerModel::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx)
{
}

void BaseEMTPlantControllerModel::init(gridpack::RealType *values)
{
}

bool BaseEMTPlantControllerModel::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

void BaseEMTPlantControllerModel::write(const char* signal, char* string)
{
}

/**
 * Set Event
 */
void BaseEMTPlantControllerModel::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{

}

int BaseEMTPlantControllerModel::matrixNumValues()
{
  return 0;
}

void BaseEMTPlantControllerModel::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  *nvals = 0;
}

void BaseEMTPlantControllerModel::vectorGetValues(gridpack::RealType *values)
{
}

void BaseEMTPlantControllerModel::setValues(gridpack::RealType *values)
{
}
