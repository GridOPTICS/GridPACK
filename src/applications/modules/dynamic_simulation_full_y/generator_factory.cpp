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

#include "generator_factory.hpp"
#include "classical.hpp"
#include "gensal.hpp"
#include "genrou.hpp"
#include "wsieg1.hpp"
#include "exdc1.hpp"
#include "esst1a.hpp"
#include "wshygp.hpp"
#include "psssim.hpp"

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
  } else if (type == "GENSAL") {
    gridpack::dynamic_simulation::GensalGenerator *tmp;
    tmp =  new gridpack::dynamic_simulation::GensalGenerator;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseGeneratorModel*>(tmp);
  } else if (type == "GENROU") {
    gridpack::dynamic_simulation::GenrouGenerator *tmp;
    tmp =  new gridpack::dynamic_simulation::GenrouGenerator;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseGeneratorModel*>(tmp);
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
  if (type == "EXDC1") {
    gridpack::dynamic_simulation::Exdc1Model *tmp;
    tmp =  new gridpack::dynamic_simulation::Exdc1Model;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseExciterModel*>(tmp);
  } else if (type == "ESST1A") {
    gridpack::dynamic_simulation::Esst1aModel *tmp;
    tmp =  new gridpack::dynamic_simulation::Esst1aModel;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseExciterModel*>(tmp);
  } else {
    ret = NULL;
  }
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
  if (type == "WSIEG1") {
    gridpack::dynamic_simulation::Wsieg1Model *tmp;
    tmp =  new gridpack::dynamic_simulation::Wsieg1Model;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseGovernorModel*>(tmp);
  } else if (type == "WSHYGP") {
    gridpack::dynamic_simulation::WshygpModel *tmp;
    tmp =  new gridpack::dynamic_simulation::WshygpModel;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseGovernorModel*>(tmp);
  } else {
    ret = NULL;
  }
  return ret;

}

/**
 * Create a PSS model
 * @param model string containing model type
 * @return pointer to model. If string does not correspond to a PSS
 * model, then return NULL pointer
 */
gridpack::dynamic_simulation::BasePssModel*
gridpack::dynamic_simulation::GeneratorFactory::createPssModel(
    std::string model)
{
  std::string type = p_util.trimQuotes(model);
  p_util.toUpper(type);

  gridpack::dynamic_simulation::BasePssModel* ret;
  if (type == "PSSSIM") {
    gridpack::dynamic_simulation::PsssimModel *tmp;
    tmp =  new gridpack::dynamic_simulation::PsssimModel;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BasePssModel*>(tmp);
  } else {
    ret = NULL;
  }
  return ret;

}
