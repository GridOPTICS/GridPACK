/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ca_driver.hpp
 * @author Bruce Palmer
 * @date   February 10, 2014
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _ca_driver_h_
#define _ca_driver_h_

#include "gridpack/configuration/configuration.hpp"

namespace gridpack {
namespace contingency_analysis {

struct Contingency
{
  std::string p_name;
  int p_id;
  std::vector<int> p_from;
  std::vector<int> p_to;
  std::vector<std::string> p_ckt;
};

// Calling program for contingency analysis application
class CADriver
{
  public:
    /**
     * Basic constructor
     */
    CADriver(void);

    /**
     * Basic destructor
     */
    ~CADriver(void);

    /**
     * Get list of contingencies from external file
     * @param cursor pointer to contingencies in input deck
     * @return vector of contingencies
     */
    std::vector<gridpack::contingency_analysis::Contingency> getContingencies(
        gridpack::utility::Configuration::CursorPtr &cursor);


    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

    private:
};

} // contingency analysis 
} // gridpack
#endif
