// -------------------------------------------------------------
/**
 * @file   pf_factory.hpp
 * @author Bruce Palmer
 * @date   July 1, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _dynsim_factory_h_
#define _dynsim_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/factory/base_factory.hpp"

namespace gridpack {
namespace dynsim {

class PFFactory
  : gridpack::factory::BaseFactory {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    PFFactory(boost::shared_ptr<gridpack::network::BaseNetwork<gridpack::component::BaseBusComponent,
                        gridpack::component::BaseBranchComponent> > network);

    /**
     * Basic destructor
     */
    ~PFFactory();

    /**
     * Create the admittance (Y-Bus) matrix
     */
    void setYBus(void);

  private:

    boost::shared_ptr<gridpack::network::BaseNetwork<gridpack::component::BaseBusComponent,
                              gridpack::component::BaseBranchComponent> > p_network;
};

} // dynsim
} // gridpack
#endif
