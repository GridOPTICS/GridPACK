/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   load_factory.cpp
 * @author Renke Huang
 * @Last modified:   August 01, 2016
 * 
 * @brief  
 * 
 * 
 */

#include "load_factory.hpp"
#include "acmotor.hpp"
#include "motorw.hpp"
#include "ieel.hpp"


/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::LoadFactory::LoadFactory(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::LoadFactory::~LoadFactory(void)
{
}

/**
 * Create a load model
 * @param model string containing model type
 * @return pointer to model. If string does not correspond
 * to a load model, then return NULL pointer
 */
gridpack::dynamic_simulation::BaseLoadModel*
  gridpack::dynamic_simulation::LoadFactory::createLoadModel(
    std::string model)
{
  std::string type = p_util.trimQuotes(model);
  p_util.toUpper(type);

  gridpack::dynamic_simulation::BaseLoadModel* ret;
  if (type == "ACMTBLU1") { //SJIN: check definition of types
    gridpack::dynamic_simulation::AcmotorLoad *tmp;
    tmp =  new gridpack::dynamic_simulation::AcmotorLoad;
    ret = dynamic_cast<gridpack::dynamic_simulation::BaseLoadModel*>(tmp);
  } else if (type == "IEEL") {
    gridpack::dynamic_simulation::IeelLoad *tmp;
    tmp =  new gridpack::dynamic_simulation::IeelLoad;
    ret = dynamic_cast<gridpack::dynamic_simulation::BaseLoadModel*>(tmp);
  } else if (type == "MOTORW") {
	gridpack::dynamic_simulation::MotorwLoad *tmp;
    tmp =  new gridpack::dynamic_simulation::MotorwLoad;
    ret = dynamic_cast<gridpack::dynamic_simulation::BaseLoadModel*>(tmp);
  } else if (type == "CIM6BL") {
	gridpack::dynamic_simulation::MotorwLoad *tmp;
    tmp =  new gridpack::dynamic_simulation::MotorwLoad;
    ret = dynamic_cast<gridpack::dynamic_simulation::BaseLoadModel*>(tmp);
  } else {
    ret = NULL;
  }
  return ret;
}


