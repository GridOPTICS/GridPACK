// -------------------------------------------------------------
/**
 * @file   pf_app.cpp
 * @author Bruce Palmer
 * @date   2013-09-09 15:17:45 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/math/linear_solver.hpp"
#include "gridpack/math/nonlinear_solver.hpp"
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
  gridpack::utility::Configuration config;
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
  //factory.setGBus();

  // make Sbus components to create S vector
  factory.setSBus();

  /*gridpack::mapper::BusVectorMap<PFNetwork> vMap(network);
  boost::shared_ptr<gridpack::math::Vector> S = vMap.mapToVector();
  S->print();*/

  // Set PQ
  //factory.setPQ();
  factory.setMode(RHS); 
  gridpack::mapper::BusVectorMap<PFNetwork> vMap(network);
  boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
  PQ->print();

  factory.setMode(Jacobian);
  gridpack::mapper::FullMatrixMap<PFNetwork> jMap(network);
  boost::shared_ptr<gridpack::math::Matrix> J = jMap.mapToMatrix();
  J->print(); 

  //factory.setMode(SetX);
  //boost::shared_ptr<gridpack::math::Vector> X(PQ->clone());
  //X->print(); 

  // Set up bus data exchange buffers. Need to decide what data needs to be
  // exchanged
  factory.setExchange();

  // Create bus data exchange
  network->initBusUpdate();
  network->updateBuses();

  // FIXME: how does one obtain the current network state (voltage/phase in this case)?
  // get the current network state vector
  // gridpack::mapper::BusVectorMap<PFNetwork> vMap(network);
  // factory.setState(); // or something

#if 0
  // Need to make sure that there is a mode for creating the X vector
  factory.setMode(RHS); 
  gridpack::mapper::BusVectorMap<PFNetwork> xMap(network);
  boost::shared_ptr<gridpack::math::Vector> X = xMap.mapToVector();


  gridpack::math::FunctionBuilder fbuilder = factory;
  gridpack::math::JacobianBuilder jbuilder = factory;
  gridpack::math::NonlinearSolver solver(network->communicator(),
                                         X->local_size(), jbuilder, fbuilder);
  
  solver.solve(*X);

#endif

#if 0
  // Initial PQ matrix
  factory.setPQ();
  gridpack::mapper::BusVectorMap<PFNetwork> vMap(network);
  boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
  PQ->print();

  boost::shared_ptr<gridpack::math::Vector> X(PQ->clone());

  gridpack::mapper::FullMatrixMap<PFNetwork> jMap(network);
  boost::shared_ptr<gridpack::math::Matrix> J = jMap.mapToMatrix();

  // Set up bus data exchange buffers. Need to decide what data needs to be
  // exchanged
  factory.setExchange();

  // Create bus data exchange
  network->initBusUpdate();
  network->updateBuses();

  // Convergence and iteration parameters
  double tolerance;
  int max_iteration;
  ComplexType tol;

  // These need to eventually be set using configuration file
  tolerance = 1.0e-6;
  max_iteration = 100;

  // Create linear solver
  gridpack::math::LinearSolver isolver(*J);

  tol = 2.0*tolerance;
  int iter = 0;

  // First iteration
  X->zero(); //might not need to do this
  isolver.solve(*PQ, *X);
  tol = X->norm2();

  while (real(tol) <= tolerance && iter <=max_iteration) {
    // Push current values in X vector back into network components
    // Need to implement setValues method in PFBus class in order for this to
    // work
    vMap.mapToBus(X);

    // Exchange data between ghost buses (I don't think we need to exchange data
    // between branches)
    network->updateBuses();

    // Create new versions of Jacobian and PQ vector
    vMap.mapToVector(PQ);
    jMap.mapToMatrix(J);

    // Create linear solver
    gridpack::math::LinearSolver solver(*J);
    X->zero(); //might not need to do this
    solver.solve(*PQ, *X);

    tol = X->norm2();
    iter++;
  }
#endif

  // Set Jacobian matrix
  // Why "getJacobian method is in PFBranch?
  // Chen 8_27_2013
  /*factory.setJacobian(); 

  gridpack::mapper::FullMatrixMap<PFNetwork> mMap(network);
  boost::shared_ptr<gridpack::math::Matrix> J = mMap.mapToMatrix();
  J->print();*/

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
