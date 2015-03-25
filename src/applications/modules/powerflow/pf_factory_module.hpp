/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_factory_module.hpp
 * @author Bruce Palmer
 * @date   2014-01-28 11:33:42 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _pf_factory_module_h_
#define _pf_factory_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/components/pf_matrix/pf_components.hpp"

namespace gridpack {
namespace powerflow {

/// The type of network used in the powerflow application
typedef network::BaseNetwork<PFBus, PFBranch > PFNetwork;

class PFFactoryModule
  : public gridpack::factory::BaseFactory<PFNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    PFFactoryModule(NetworkPtr network);

    /**
     * Basic destructor
     */
    ~PFFactoryModule();

    /**
     * Create the admittance (Y-Bus) matrix
     */
    void setYBus(void);

    /**
     * Make SBus vector 
     */
    void setSBus(void);

    /**
      * Update pg of specified bus element based on their genID
      * @param busID
      * @param genID
      * @param value
      */
    void updatePg(int busID, std::string genID, double value);

    /**
     * Create the PQ 
     */
    void setPQ(void);

  private:

    NetworkPtr p_network;
};

} // powerflow
} // gridpack
#endif
