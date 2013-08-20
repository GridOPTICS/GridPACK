// -------------------------------------------------------------
/**
 * @file   pf_factory.hpp
 * @author Bruce Palmer
 * @date   2013-08-08 10:20:30 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _pf_factory_h_
#define _pf_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/applications/powerflow/pf_components.hpp"

namespace gridpack {
namespace powerflow {

class PFFactory
  : public gridpack::factory::BaseFactory<PFNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    PFFactory(NetworkPtr network);

    /**
     * Basic destructor
     */
    ~PFFactory();

    /**
     * Create the admittance (Y-Bus) matrix
     */
    void setYBus(void);

    /**
     * Create SBus matrix
     */
    void setSBus(void);

    /**
     * Create the Jacobian matrix
     */
    void setJacobian(void);

    /**
     * Create RHS
     */

    void setRHS(void);

  private:

    NetworkPtr p_network;
};

} // powerflow
} // gridpack
#endif
