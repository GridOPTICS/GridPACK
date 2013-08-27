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

//enum PFMode{YBus, Jacobian};

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
     * Load data
     */ 
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

    /**
     * Set the mode to control what matrices and vectors are built when using the mapper 
     */
    void setMode(int mode);

    /**
     * Create the admittance (Y-Bus) matrix
     */
    void setYBus(void);

    /**
     * Find GBus vector 
     */
    void setGBus(void);

    /**
     * Make SBus vector 
     */
    void setSBus(void);

    /**
     * Create the PQ 
     */
    void setPQ(void);

    /**
     * Create the Jacobian matrix
     */
    void setJacobian(void);

  private:

    NetworkPtr p_network;
};

} // powerflow
} // gridpack
#endif
