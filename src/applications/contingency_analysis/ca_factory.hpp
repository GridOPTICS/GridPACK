/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ca_factory.hpp
 * @author Yousu Chen 
 * @date   January 20, 2014
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _ca_factory_h_
#define _ca_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/applications/contingency_analysis/ca_components.hpp"
#include "gridpack/math/matrix.hpp"

namespace gridpack {
namespace contingency_analysis {

class CAFactory
  : public gridpack::factory::BaseFactory<CANetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    CAFactory(NetworkPtr network);

    /**
     * Basic destructor
     */
    ~CAFactory();

    /**
     * Create the admittance (Y-Bus) matrix
     */
    void setYBus(void);

  private:

    NetworkPtr p_network;
};

} // contingency_analysis 
} // gridpack
#endif
