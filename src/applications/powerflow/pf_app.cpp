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
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/parser/Parser.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/applications/powerflow/pf_factory.hpp"


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

  // read configuration file
  gridpack::utilities::Configuration config;
//  config.open("config.txt", world);

  // load input file
  gridpack::parser::PTI23_parser<PFNetwork> parser(network);
 // parser.getCase("118_pti_v29.raw");
  parser.getCase("IEEE14.raw");
  parser.createNetwork();

  // partition network
  network->partition();

  // create factory
  gridpack::powerflow::PFFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();

  // set YBus components so that you can create Y matrix
  factory.setYBus();

  factory.setMode(YBus); 
 
  gridpack::mapper::FullMatrixMap<PFNetwork> mMap(network);
  boost::shared_ptr<gridpack::math::Matrix> Y = mMap.mapToMatrix();
  Y->print();

  // set GBus components to create G vector 
  factory.setGBus();

  // make Sbus components to create S vector
  factory.setSBus();

  /*gridpack::mapper::BusVectorMap<PFNetwork> vMap(network);
  boost::shared_ptr<gridpack::math::Vector> S = vMap.mapToVector();
  S->print();*/

  // Set PQ
  //factory.setPQ();

  gridpack::mapper::BusVectorMap<PFNetwork> vMap(network);
  boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
  PQ->print();

  factory.setMode(Jacobian);

  // Set Jacobian matrix
  // Why "getJacobian method is in PFBranch?
  // Chen 8_27_2013
//  factory.setJacobian(); 

/*  // Start AC N-R Solver

  // FIND the first mismatch
  // MIS: vector
  // V: vector
  // YBus: Matrix
  // SBUS: vector
  // MIS = V * conj (YBus * V) - SBUS
  factory.calMis();

  // Assume that matrix A and vector V have been properly set up and start
  // creating solver loop.
  factory.setMode(YBus);
  gridpack::mapper::FullMatrixMap<PFNetwork> mMap(network);
  boost::shared_ptr<gridpack::math::Matrix> Y = mMap.mapToMatrix();
  Y->print();
  return;

  gridpack::mapper::BusVectorMap<PFNetwork> vMap(network);
  boost::shared_ptr<gridpack::math::Vector> V = vMap.mapToVector();

  boost::shared_ptr<gridpack::math::Vector> SBus(V->clone());
  boost::shared_ptr<gridpack::math::Vector> MIS(V->clone());
*/
}
