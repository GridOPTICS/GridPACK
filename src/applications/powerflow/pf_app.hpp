// -------------------------------------------------------------
/**
 * @file   pf_app.hpp
 * @author Bruce Palmer
 * @date   July 23, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _pf_app_h_
#define _pf_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/applications/powerflow/pf_components.hpp"
#include "gridpack/applications/powerflow/pf_factory.hpp"

typedef gridpack::network::BaseNetwork<
        gridpack::powerflow::PFBus,
        gridpack::powerflow::PFBranch> PFNetwork;

namespace gridpack {
namespace powerflow {

// Calling program for powerflow application

class PFApp
{
  public:
    /**
     * Basic constructor
     */
    PFApp(void);

    /**
     * Basic destructor
     */
    ~PFApp(void);

    /**
     * Execute application
     */
    void execute(void);

  private:
};

} // powerflow
} // gridpack
#endif
