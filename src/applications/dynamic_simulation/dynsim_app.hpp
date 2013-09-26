// -------------------------------------------------------------
/**
 * @file   dynsim_app.hpp
 * @author Shuangshuang Jin
 * @date   September 19, 2013
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _dynsim_app_h_
#define _dynsim_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/applications/dynamis_simulation/dynsim_factory.hpp"

namespace gridpack {
namespace dynamic_simulation {

// Calling program for dynamic simulation application

class DynSimApp
{
  public::
    /**
     * Basic constructor
     */
     DynSimApp(void);

    /**
     * Basic destructor
     */
    ~DynSimApp(void);

    /**
     * Execute application
     */
    void execute(void);
 
    private:
};

} // dynamic simulation
} // gridpack
#endif
