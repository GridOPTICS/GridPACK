/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsf_generator_factory_module.hpp
 * @author Bruce Palmer
 * @Last modified:   May 19, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef dsf_generator_factory_module_h_
#define dsf_generator_factory_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/include/gridpack.hpp"
#include "base_classes/base_generator_model.hpp"
#include "base_classes/base_exciter_model.hpp"
#include "base_classes/base_governor_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class DSFGeneratorFactory
{
  public:
    /**
     * Basic constructor
     */
    DSFGeneratorFactory();

    /**
     * Basic destructor
     */
    virtual ~DSFGeneratorFactory();

    /**
     * Create a generator model
     * @param model string containing model type
     * @return pointer to model. If string does not correspond to a generator
     * model, then return NULL pointer
     */
    DSFBaseGeneratorModel* createGeneratorModel(std::string model);

    /**
     * Create an exciter model
     * @param model string containing model type
     * @return pointer to model. If string does not correspond to a exciter
     * model, then return NULL pointer
     */
    DSFBaseExciterModel* createExciterModel(std::string model);

    /**
     * Create an governor model
     * @param model string containing model type
     * @return pointer to model. If string does not correspond to a governor
     * model, then return NULL pointer
     */
    DSFBaseGovernorModel* createGovernorModel(std::string model);

  private:

    gridpack::utility::StringUtils p_util;
};
}  // dynamic_simulation
}  // gridpack
#endif
