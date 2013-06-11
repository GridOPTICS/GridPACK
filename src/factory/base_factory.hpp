// -------------------------------------------------------------
/**
 * @file   base_factory.hpp
 * @author Bruce Palmer
 * @date   June 3, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _base_factory_h_
#define _base_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"

// Base factory class that contains functions that are generic to all
// applications.

namespace gridpack{
namespace factory{

class BaseFactory {
  public:
    /**
     * Constructor
     */
    BaseFactory(void);

    /**
     * Destructor
     */
    virtual ~BaseFactory(void);

    /**
     * Set pointers in each bus and branch component so that it points to
     * connected buses and branches. This routine operates on the generic
     * BaseBusComponent and BaseBranchComponent interfaces.
     * @param network: The network that contains the components that need to be
     * set
     */
    virtual void setComponents(gridpack::network::BaseNetwork<gridpack::component::BaseBusComponent>,
                               gridpack::component::BaseBranchComponent> >
                               *network);

  private:

};

}    // factory
}    // gridpack
#endif
