/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
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
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/applications/dynamic_simulation/ds_factory.hpp"

namespace gridpack {
namespace dynamic_simulation {

// Small utility structure to encapsulate information about fault events
struct Event{
  double start,end;
  int from_idx, to_idx;
  double step;
};


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
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

    /**
     * Utility function to convert faults that are in event list into
     * internal data structure that can be used by code
     * @param cursors list of cursors pointing to individual events in input
     * deck
     * @return list of event data structures
     */
    std::vector<Event>
      setFaultEvents(std::vector<gridpack::utility::Configuration::CursorPtr > cursors);
 
    private:
};

} // dynamic simulation
} // gridpack
#endif
