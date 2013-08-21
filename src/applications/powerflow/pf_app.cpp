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
#include "gridpack/mapper/bus_vector_map.hpp"
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
  config.open("config.txt", world);

  // load input file
  gridpack::parser::PTI23_parser parser;
  parser.getCase("PTI23file");

  // partition network
  network->partition();

  // create factory
  gridpack::powerflow::PFFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();

  // set YBus components so that you can create Y matrix
  factory.setYBus();
#if 0
  // SetYbus
  gridpack::powerflow::pf_factory factory;
  factory.setGBus();
  factory.setSBus();
  factory.setInit();

  // Start AC N-R Solver

  // FIND the first mismatch
  // MIS: vector
  // V: vector
  // YBus: Matrix
  // SBUS: vector
  // MIS = V * conj (YBus * V) - SBUS
  factory.calMis();
#endif

  // Assume that matrix A and vector V have been properly set up and start
  // creating solver loop.
  factory.setMode(YBus);
  gridpack::mapper::FullMatrixMap<PFNetwork> mMap(network);
  boost::shared_ptr<gridpack::math::Matrix> Y = mMap.mapToMatrix();
  Y->print();

  gridpack::mapper::BusVectorMap<PFNetwork> vMap(network);
  boost::shared_ptr<gridpack::math::Vector> V = vMap.mapToVector();
}
