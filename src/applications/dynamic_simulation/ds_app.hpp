// -------------------------------------------------------------
/**
 * @file   ds_app.hpp
 * @author Shuangshuang Jin
 * @date   September 19, 2013
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _ds_app_h_
#define _ds_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/applications/dynamic_simulation/ds_factory.hpp"

namespace gridpack {
namespace dynamic_simulation {

// Calling program for dynamic simulation application

class DSApp
{
  public:
    /**
     * Basic constructor
     */
     DSApp(void);

    /**
     * Basic destructor
     */
    ~DSApp(void);

    /**
     * Execute application
     */
    void execute(void);
 
    private:
};

} // dynamic simulation
} // gridpack
#endif
