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

#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/applications/powerflow/pf_app.hpp"
#include "gridpack/mapper/full_map.hpp"

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

  // Assume that matrix A and vector B have been properly set up and start
  // creating solver loop.
  gridpack::mapper::FullMatrixMap<PFNetwork> map(network);
  boost::shared_ptr<gridpack::math::Matrix> A = map.mapToMatrix();
  // boost::shared_ptr<gridpack::math::Vector> B(new gridpack::math::Vector);
  // boost::shared_ptr<gridpack::math::Vector> X(new gridpack::math::Vector);
}
