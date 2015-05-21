/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   uc_app.hpp
 * @author 
 * @date   
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _uc_app_h_
#define _uc_app_h_

namespace gridpack {
namespace unit_commitment {

// Calling program for unit commitment application

class UCApp
{
  public:
    /**
     * Basic constructor
     */
    UCApp(void);

    /**
     * Basic destructor
     */
    ~UCApp(void);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);


  private:
};

} // unit_commitment
} // gridpack
#endif
