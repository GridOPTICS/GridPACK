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
#include "gdform.hpp"
#include "regca1.hpp"
#include "regcb1.hpp"
#include "regcc1.hpp"
#include "reeca1.hpp"
#include "repca1.hpp"
#include "wsieg1.hpp"
#include "exdc1.hpp"
#include "ieeet1.hpp"
#include "esst1a.hpp"
#include "wshygp.hpp"
#include "psssim.hpp"
#include "tgov1.hpp"
#include "sexs.hpp"
#include "gast.hpp"
#include "hygov.hpp"
#include "wtara1.hpp"
#include "wtdta1.hpp"
#include "wtpta1.hpp"
#include "wttqa1.hpp"
#include "stab1.hpp"

#include <stdio.h>

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
  } else if (type == "GDFORM") {
    gridpack::dynamic_simulation::GridFormingGenerator *tmp;
    tmp =  new gridpack::dynamic_simulation::GridFormingGenerator;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseGeneratorModel*>(tmp);
  } else if (type == "REGCA1") {
    gridpack::dynamic_simulation::Regca1Generator *tmp;
    tmp =  new gridpack::dynamic_simulation::Regca1Generator;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseGeneratorModel*>(tmp);
  } else if (type == "REGCB1") {
    gridpack::dynamic_simulation::Regcb1Generator *tmp;
    tmp =  new gridpack::dynamic_simulation::Regcb1Generator;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseGeneratorModel*>(tmp);
  } else if (type == "REGCC1") {
    gridpack::dynamic_simulation::Regcc1Generator *tmp;
    tmp =  new gridpack::dynamic_simulation::Regcc1Generator;
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
  }  else if (type == "IEEET1") {
    gridpack::dynamic_simulation::Ieeet1Model *tmp;
    tmp =  new gridpack::dynamic_simulation::Ieeet1Model;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseExciterModel*>(tmp);
  } else if (type == "SEXS") {
      gridpack::dynamic_simulation::SexsModel *tmp;
      tmp =  new gridpack::dynamic_simulation::SexsModel;
      ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseExciterModel*>(tmp);
  } else if (type == "ESST1A") {
    gridpack::dynamic_simulation::Esst1aModel *tmp;
    tmp =  new gridpack::dynamic_simulation::Esst1aModel;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseExciterModel*>(tmp);
  }else if (type == "REECA1") {
    gridpack::dynamic_simulation::Reeca1Model *tmp;
    tmp =  new gridpack::dynamic_simulation::Reeca1Model;
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
  } else if (type == "TGOV1") {
    gridpack::dynamic_simulation::Tgov1Model *tmp;
    tmp =  new gridpack::dynamic_simulation::Tgov1Model;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseGovernorModel*>(tmp);
  } else if (type == "GAST") {
    gridpack::dynamic_simulation::GastModel *tmp;
    tmp =  new gridpack::dynamic_simulation::GastModel;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseGovernorModel*>(tmp);
  } else if (type == "HYGOV") {
    gridpack::dynamic_simulation::HygovModel *tmp;
    tmp =  new gridpack::dynamic_simulation::HygovModel;
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
  } else if (type == "STAB1") {
    gridpack::dynamic_simulation::Stab1Model *tmp;
    tmp =  new gridpack::dynamic_simulation::Stab1Model;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BasePssModel*>(tmp);
  } else {
    ret = NULL;
  }
  return ret;

}


/**
 * Create a plant controller model
 * @param model string containing model type
 * @return pointer to model. If string does not correspond to a PSS
 * model, then return NULL pointer
 */
gridpack::dynamic_simulation::BasePlantControllerModel*
gridpack::dynamic_simulation::GeneratorFactory::createPlantControllerModel(
    std::string model)
{
  std::string type = p_util.trimQuotes(model);
  p_util.toUpper(type);

  gridpack::dynamic_simulation::BasePlantControllerModel* ret;
  if (type == "REPCA1" || type == "REPCTA1") {
    gridpack::dynamic_simulation::Repca1Model *tmp;
    tmp =  new gridpack::dynamic_simulation::Repca1Model;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BasePlantControllerModel*>(tmp);
  } else {
    ret = NULL;
  }
  return ret;

}

/**
 * Create a wind mechanical model
 * @param model string containing model type
 * @return pointer to model. If string does not correspond to a wind mechanical model
 * model, then return NULL pointer
 */
gridpack::dynamic_simulation::BaseMechanicalModel*
gridpack::dynamic_simulation::GeneratorFactory::createMechanicalModel(
    std::string model)
{
  std::string type = p_util.trimQuotes(model);
  p_util.toUpper(type);

  gridpack::dynamic_simulation::BaseMechanicalModel* ret;
  if (type == "WTARA1" ) {
    gridpack::dynamic_simulation::Wtara1Model *tmp;
    tmp =  new gridpack::dynamic_simulation::Wtara1Model;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseMechanicalModel*>(tmp);
  } else if(type == "WTDTA1" ) {
    gridpack::dynamic_simulation::Wtdta1Model *tmp;
    tmp =  new gridpack::dynamic_simulation::Wtdta1Model;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseMechanicalModel*>(tmp);

  } else if(type == "WTPTA1" ) {
    gridpack::dynamic_simulation::Wtpta1Model *tmp;
    tmp =  new gridpack::dynamic_simulation::Wtpta1Model;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseMechanicalModel*>(tmp);

  } else if(type == "WTTQA1" ) {
    gridpack::dynamic_simulation::Wttqa1Model *tmp;
    tmp =  new gridpack::dynamic_simulation::Wttqa1Model;
    ret =
      dynamic_cast<gridpack::dynamic_simulation::BaseMechanicalModel*>(tmp);
    
  } else  ret = NULL;
  
  return ret;

}

