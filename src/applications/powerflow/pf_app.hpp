// -------------------------------------------------------------
/**
 * @file   pf_app.hpp
 * @author Bruce Palmer
 * @date   2013-08-08 10:18:32 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _pf_app_h_
#define _pf_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/applications/powerflow/pf_factory.hpp"

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
