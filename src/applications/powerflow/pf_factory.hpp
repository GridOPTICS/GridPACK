// -------------------------------------------------------------
/**
 * @file   pf_factory.hpp
 * @author Bruce Palmer
 * @date   2013-09-09 14:25:24 d3g096
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
#include "gridpack/math/matrix.hpp"

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

    /// Operator to make this compatible with math::JacobianBuilder
  void operator() (const math::Vector& x, math::Matrix& J);

    /// Operator to make this compatible with math::FunctionBuilder
  void operator() (const math::Vector& x, math::Vector& F);

  private:

    NetworkPtr p_network;
};

} // powerflow
} // gridpack
#endif
