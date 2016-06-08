/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   dsf_generator_factory_module.cpp
 * @author Bruce Palmer
 * @Last modified:   May 18, 2015
 * 
 * @brief  
 * 
 * 
 */

#include "gridpack/include/gridpack.hpp"
#include "dsf_generator_factory_module.hpp"
#include "model_classes/classical.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::DSFGeneratorFactory::DSFGeneratorFactory(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::DSFGeneratorFactory::~DSFGeneratorFactory(void)
{
}

/**
 * Create a generator model
 * @param model string containing model type
 * @return pointer to model. If string does not correspond
 * to a generator model, then return NULL pointer
 */
gridpack::dynamic_simulation::DSFBaseGeneratorModel*
  gridpack::dynamic_simulation::DSFGeneratorFactory::createGeneratorModel(
    std::string model)
{
  std::string type = p_util.trimQuotes(model);
  p_util.toUpper(type);

  gridpack::dynamic_simulation::DSFBaseGeneratorModel* ret;
  if (type == "GENCLS") {
    gridpack::dynamic_simulation::DSFClassicalGenerator *tmp;
    tmp =  new gridpack::dynamic_simulation::DSFClassicalGenerator;
    ret = dynamic_cast<gridpack::dynamic_simulation::DSFBaseGeneratorModel*>(tmp);
  } else {
    ret = NULL;
  }
  return ret;
}

/**
 * Create an exciter model
 * @param model string containing model type
 * @return pointer to model. If string does not correspond to a exciter
 * model, then return NULL pointer
 */
gridpack::dynamic_simulation::DSFBaseExciterModel*
gridpack::dynamic_simulation::DSFGeneratorFactory::createExciterModel(
    std::string model)
{
  std::string type = p_util.trimQuotes(model);
  p_util.toUpper(type);

  gridpack::dynamic_simulation::DSFBaseExciterModel* ret;
  ret = NULL;
  return ret;
}

/**
 * Create a governor model
 * @param model string containing model type
 * @return pointer to model. If string does not correspond to a governor
 * model, then return NULL pointer
 */
gridpack::dynamic_simulation::DSFBaseGovernorModel*
gridpack::dynamic_simulation::DSFGeneratorFactory::createGovernorModel(
    std::string model)
{
  std::string type = p_util.trimQuotes(model);
  p_util.toUpper(type);

  gridpack::dynamic_simulation::DSFBaseGovernorModel* ret;
  ret = NULL;
  return ret;
}
