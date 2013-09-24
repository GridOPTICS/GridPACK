// -------------------------------------------------------------
/**
 * @file   dynsim_factory.hpp
 * @author Shuangshuang Jin 
 * @date   September 19, 2013
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
#include "gridpack/applications/dynamic_simulation/dynsim_components.hpp"
#include "gridpack/math/matrix.hpp"

namespace gridpack {
namespace dynsim {

class DynSimFactory
  : public gridpack::factory::BaseFactory<DynSimNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    DynSimFactory(NetworkPtr network);

    /**
     * Basic destructor
     */
    ~DynSimFactory();

    /**
     * Create the admittance (Y-Bus) matrix
     */
    void setYBus(void);

  private:

    NetworkPtr p_network;
};

} // dynsim
} // gridpack
#endif
