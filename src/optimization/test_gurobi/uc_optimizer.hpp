/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   uc_optimizer.hpp
 * @author 
 * @date   
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _uc_optimizer_h_
#define _uc_optimizer_h_

#include "gridpack/optimization/optimization.hpp"
#include "boost/smart_ptr/shared_ptr.hpp"
#include "uc_components.hpp"

namespace gridpack {
namespace unit_commitment {

// This example only needs the functionality in the base optimization class

class UCoptimizer
    : public gridpack::optimization::Optimizer<UCNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with optimization
     */
    UCoptimizer(boost::shared_ptr<UCNetwork> network) 
      : gridpack::optimization::Optimizer<UCNetwork>(network)
    {
      p_network = network;
    }


    /**
     * Basic destructor
     */
    ~UCoptimizer() {}

  private:

    NetworkPtr p_network;
};

} // uc_commitment
} // gridpack
#endif
