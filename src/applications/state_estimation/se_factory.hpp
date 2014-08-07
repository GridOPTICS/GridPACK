/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_factory.hpp
 * @author Yousu Chen 
 * @date   2/24/2014 
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _se_factory_h_
#define _se_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/include/gridpack.hpp"
#include "se_components.hpp"

namespace gridpack {
namespace state_estimation {

class SEFactory
  : public gridpack::factory::BaseFactory<SENetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    SEFactory(NetworkPtr network);

    /**
     * Basic destructor
     */
    ~SEFactory();

    /**
     * Create the admittance (Y-Bus) matrix
     */
    void setYBus(void);

    /**
     * Disribute measurements
     * @param measurements a vector containing all measurements
     */
    void setMeasurements(std::vector<Measurement> measurements);

  private:

    NetworkPtr p_network;
};

} // state_estimation
} // gridpack
#endif
