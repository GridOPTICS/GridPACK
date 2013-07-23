// -------------------------------------------------------------
/**
 * @file   pf_main.cpp
 * @author Bruce Palmer
 * @date   July 23, 2013
 * 
 * @brief  
 */
// -------------------------------------------------------------

#include "gridpack/applications/powerflow/pf_app.hpp"

// Calling program for the powerflow applications

main(int *argc, char **argv)
{
  gridpack::powerflow::PFApp app;
  app.execute();
}
