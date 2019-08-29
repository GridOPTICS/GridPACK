/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   load_factory.hpp
 * @author Renke Huang
 * @Last modified:   August 01, 2016
 * 
 * @brief  
 * 
 * 
 */

#ifndef load_factory_h_
#define load_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_classes/base_load_model.hpp"
#include "gridpack/utilities/string_utils.hpp"
//#include "base_exciter_model.hpp"
//#include "base_governor_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class LoadFactory
{
  public:
    /**
     * Basic constructor
     */
    LoadFactory();

    /**
     * Basic destructor
     */
    virtual ~LoadFactory();

    /**
     * Create a load model
     * @param model string containing model type
     * @return pointer to model. If string does not correspond to a load
     * model, then return NULL pointer
     */
    BaseLoadModel* createLoadModel(std::string model);

  private:

    gridpack::utility::StringUtils p_util;
};
}  // dynamic_simulation
}  // gridpack
#endif
