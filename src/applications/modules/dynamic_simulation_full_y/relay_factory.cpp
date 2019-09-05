/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   relay_factory.cpp
 * @author Renke Huang
 * @Last modified:   August 01, 2016
 * 
 * @brief  
 * 
 * 
 */

#include "relay_factory.hpp"
#include "lvshbl.hpp"
#include "frqtpat.hpp"
#include "distr1.hpp"


/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::RelayFactory::RelayFactory(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::RelayFactory::~RelayFactory(void)
{
}

/**
 * Create a relay model
 * @param model string containing model type
 * @return pointer to model. If string does not correspond
 * to a relay model, then return NULL pointer
 */
gridpack::dynamic_simulation::BaseRelayModel*
  gridpack::dynamic_simulation::RelayFactory::createRelayModel(
    std::string model)
{
  std::string type = p_util.trimQuotes(model);
  p_util.toUpper(type);

  gridpack::dynamic_simulation::BaseRelayModel* ret;
  if (type == "LVSHBL") {
    gridpack::dynamic_simulation::LvshblRelay *tmp;
    tmp =  new gridpack::dynamic_simulation::LvshblRelay;
    ret = dynamic_cast<gridpack::dynamic_simulation::BaseRelayModel*>(tmp);
  } else if (type == "FRQTPAT") {
    gridpack::dynamic_simulation::FrqtpatRelay *tmp;
    tmp =  new gridpack::dynamic_simulation::FrqtpatRelay;
    ret = dynamic_cast<gridpack::dynamic_simulation::BaseRelayModel*>(tmp);
  } else if (type == "DISTR1") {
	gridpack::dynamic_simulation::Distr1Relay *tmp;
    tmp =  new gridpack::dynamic_simulation::Distr1Relay;
    ret = dynamic_cast<gridpack::dynamic_simulation::BaseRelayModel*>(tmp);
  } else {
    ret = NULL;
  }
  return ret;
}


