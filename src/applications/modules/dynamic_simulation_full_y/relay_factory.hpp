/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   relay_factory.hpp
 * @author Renke Huang
 * @Last modified:   August 01, 2016
 * 
 * @brief  
 * 
 * 
 */

#ifndef relay_factory_h_
#define relay_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/string_utils.hpp"
#include "base_classes/base_relay_model.hpp"
//#include "base_exciter_model.hpp"
//#include "base_governor_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class RelayFactory
{
  public:
    /**
     * Basic constructor
     */
    RelayFactory();

    /**
     * Basic destructor
     */
    virtual ~RelayFactory();

    /**
     * Create a relay model
     * @param model string containing model type
     * @return pointer to model. If string does not correspond to a relay
     * model, then return NULL pointer
     */
    BaseRelayModel* createRelayModel(std::string model);

  private:

    gridpack::utility::StringUtils p_util;
};
}  // dynamic_simulation
}  // gridpack
#endif
