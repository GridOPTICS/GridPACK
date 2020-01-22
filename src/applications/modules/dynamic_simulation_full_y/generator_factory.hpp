/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   generator_factory.hpp
 * @author Bruce Palmer
 * @Last modified:   May 19, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef generator_factory_h_
#define generator_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_classes/base_generator_model.hpp"
#include "base_classes/base_exciter_model.hpp"
#include "base_classes/base_governor_model.hpp"
#include "base_classes/base_pss_model.hpp"
#include "gridpack/utilities/string_utils.hpp"

namespace gridpack {
namespace dynamic_simulation {
class GeneratorFactory
{
  public:
    /**
     * Basic constructor
     */
    GeneratorFactory();

    /**
     * Basic destructor
     */
    virtual ~GeneratorFactory();

    /**
     * Create a generator model
     * @param model string containing model type
     * @return pointer to model. If string does not correspond to a generator
     * model, then return NULL pointer
     */
    BaseGeneratorModel* createGeneratorModel(std::string model);

    /**
     * Create an exciter model
     * @param model string containing model type
     * @return pointer to model. If string does not correspond to a exciter
     * model, then return NULL pointer
     */
    BaseExciterModel* createExciterModel(std::string model);

    /**
     * Create an governor model
     * @param model string containing model type
     * @return pointer to model. If string does not correspond to a governor
     * model, then return NULL pointer
     */
    BaseGovernorModel* createGovernorModel(std::string model);
	
	/**
     * Create an PSS model
     * @param model string containing model type
     * @return pointer to model. If string does not correspond to a pss
     * model, then return NULL pointer
     */
    BasePssModel* createPssModel(std::string model);

  private:

    gridpack::utility::StringUtils p_util;
};
}  // dynamic_simulation
}  // gridpack
#endif
