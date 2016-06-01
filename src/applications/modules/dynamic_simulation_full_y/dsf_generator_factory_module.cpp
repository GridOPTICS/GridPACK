/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   generator_factory.cpp
 * @author Bruce Palmer
 * @Last modified:   May 18, 2015
 * 
 * @brief  
 * 
 * 
 */

#include "gridpack/include/gridpack.hpp"
#include "dsf_generator_factory_module.hpp"
#include "classical.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::GeneratorFactory::GeneratorFactory(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::GeneratorFactory::~GeneratorFactory(void)
{
}

/**
 * Create a generator model
 * @param model string containing model type
 * @return pointer to model. If string does not correspond
 * to a generator model, then return NULL pointer
 */
gridpack::dynamic_simulation::BaseGeneratorModel*
  gridpack::dynamic_simulation::GeneratorFactory::createGeneratorModel(
    std::string model)
{
  std::string type = p_util.trimQuotes(model);
  p_util.toUpper(type);

  gridpack::dynamic_simulation::BaseGeneratorModel* ret;
  if (type == "GENCLS") {
    gridpack::dynamic_simulation::ClassicalGenerator *tmp;
    tmp =  new gridpack::dynamic_simulation::ClassicalGenerator;
    ret = dynamic_cast<gridpack::dynamic_simulation::BaseGeneratorModel*>(tmp);
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
gridpack::dynamic_simulation::BaseExciterModel*
gridpack::dynamic_simulation::GeneratorFactory::createExciterModel(
    std::string model)
{
  std::string type = p_util.trimQuotes(model);
  p_util.toUpper(type);

  gridpack::dynamic_simulation::BaseExciterModel* ret;
  ret = NULL;
  return ret;
}

/**
 * Create a governor model
 * @param model string containing model type
 * @return pointer to model. If string does not correspond to a governor
 * model, then return NULL pointer
 */
gridpack::dynamic_simulation::BaseGovernorModel*
gridpack::dynamic_simulation::GeneratorFactory::createGovernorModel(
    std::string model)
{
  std::string type = p_util.trimQuotes(model);
  p_util.toUpper(type);

  gridpack::dynamic_simulation::BaseGovernorModel* ret;
  ret = NULL;
  return ret;
}
