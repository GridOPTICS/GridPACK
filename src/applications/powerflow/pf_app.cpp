// -------------------------------------------------------------
/**
 * @file   pf_app.cpp
 * @author Bruce Palmer
 * @date   July 23, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/applications/powerflow/pf_app.hpp"

// Calling program for powerflow application

/**
 * Basic constructor
 */
gridpack::powerflow::PFApp::PFApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::powerflow::PFApp::~PFApp(void)
{
}

/**
 * Execute application
 */
void gridpack::powerflow::PFApp::execute(void)
{
  gridpack::parallel::Communicator world;
  boost::shared_ptr<PFNetwork> network(new PFNetwork(world));
}
