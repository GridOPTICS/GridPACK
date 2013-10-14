/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_app.cpp
 * @author Bruce Palmer
 * @date   2013-10-08 08:15:14 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/math/linear_solver.hpp"
#include "gridpack/math/newton_raphson_solver.hpp"
#include "gridpack/math/nonlinear_solver.hpp"
#include "gridpack/applications/powerflow/pf_app.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/applications/powerflow/pf_factory.hpp"
#include "gridpack/timer/coarse_timer.hpp"


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
  gridpack::utility::Configuration *config = gridpack::utility::Configuration::configuration();
  config->open("input.xml",world);
  gridpack::utility::Configuration::Cursor *cursor;
  cursor = config->getCursor("Configuration.Powerflow");
  std::string filename = cursor->get("networkConfiguration",
      "No network configuration specified");

  // load input file
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_pti = timer->createCategory("PTI Parser");
  timer->start(t_pti);
  gridpack::parser::PTI23_parser<PFNetwork> parser(network);
  parser.parse(filename.c_str());
  timer->stop(t_pti);

//  std::string unpartout(cursor->get("networkUnpartitionedGraph", ""));
//  std::string partout(cursor->get("networkPartitionedGraph", ""));

//  if (!unpartout.empty()) {
//    network->writeGraph(unpartout);
//  }

  // partition network
  int t_part = timer->createCategory("Partition");
  timer->start(t_part);
  network->partition();
  timer->stop(t_part);

//  if (!partout.empty()) {
//    network->writeGraph(partout);
//  }

  // create factory
  gridpack::powerflow::PFFactory factory(network);
  int t_fload = timer->createCategory("Factory Load");
  timer->start(t_fload);
  factory.load();
  timer->stop(t_fload);

  // set network components using factory
  int t_fset = timer->createCategory("Factory Set Components");
  timer->start(t_fset);
  factory.setComponents();
  timer->stop(t_fset);

  // Set up bus data exchange buffers. Need to decide what data needs to be
  // exchanged
  int t_fex = timer->createCategory("Factory Set Exchange");
  timer->start(t_fex);
  factory.setExchange();
  timer->stop(t_fex);

  // Create bus data exchange
  int t_setupdt = timer->createCategory("Set Bus Update");
  timer->start(t_setupdt);
  network->initBusUpdate();
  timer->stop(t_setupdt);


  // set YBus components so that you can create Y matrix
  factory.setYBus();

  // Create serial IO object to export data from buses
  gridpack::serial_io::SerialBusIO<PFNetwork> busIO(128,network);
  char ioBuf[128];

  int t_ymap = timer->createCategory("Create Y-matrix Mapper");
  timer->start(t_ymap);
  factory.setMode(YBus); 
  gridpack::mapper::FullMatrixMap<PFNetwork> mMap(network);
  timer->stop(t_ymap);
  int t_ybus = timer->createCategory("Create Y-matrix");
  timer->start(t_ybus);
  boost::shared_ptr<gridpack::math::Matrix> Y = mMap.mapToMatrix();
  timer->stop(t_ybus);
//  busIO.header("\nY-matrix values\n");
//  Y->print();
  Y->save("Ybus.m");

  // make Sbus components to create S vector
  factory.setSBus();
  busIO.header("\nIteration 0\n");

  // Set PQ
  int t_pqmap = timer->createCategory("Create PQ Mapper");
  timer->start(t_pqmap);
  factory.setMode(RHS); 
  gridpack::mapper::BusVectorMap<PFNetwork> vMap(network);
  timer->stop(t_pqmap);
  int t_pqvec = timer->createCategory("Create PQ Vector");
  timer->start(t_pqvec);
  boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
  timer->stop(t_pqvec);
//  busIO.header("\nPQ values\n");
//  PQ->print();
  busIO.header("\n   Elements of PQ vector\n");
  busIO.header("\n   Bus Number           P                   Q      Neighbors\n");
  busIO.write("pq");
#if 1
  factory.setMode(Jacobian);
  gridpack::mapper::FullMatrixMap<PFNetwork> jMap(network);
  boost::shared_ptr<gridpack::math::Matrix> J = jMap.mapToMatrix();
//  busIO.header("\nJacobian values\n");
//  J->print(); 

  // FIXME: how does one obtain the current network state (voltage/phase in this case)?
  // get the current network state vector
  // gridpack::mapper::BusVectorMap<PFNetwork> vMap(network);
  // factory.setState(); // or something

#if 0
  // Need to make sure that there is a mode for creating the X vector
  factory.setPQ();
  gridpack::mapper::BusVectorMap<PFNetwork> xMap(network);
  boost::shared_ptr<gridpack::math::Vector> X = xMap.mapToVector();


  gridpack::math::FunctionBuilder fbuilder = factory;
  gridpack::math::JacobianBuilder jbuilder = factory;

#if 0
  gridpack::math::NewtonRaphsonSolver solver(network->communicator(),
                                             X->local_size(), jbuilder, fbuilder);
  solver.tolerance(1.0e-06);
  solver.maximum_iterations(100);
  solver.solve(*X);
#else
  gridpack::math::NonlinearSolver solver(network->communicator(),
                                         X->local_size(), jbuilder, fbuilder);
  try {
    solver.solve(*X);
  } catch (const Exception& e) {
    std::cerr << e.what() << std::endl;
  }
#endif

  X->print();

#endif

#if 1
  boost::shared_ptr<gridpack::math::Vector> X(PQ->clone());

  // Convergence and iteration parameters
  double tolerance;
  int max_iteration;
  ComplexType tol;

  // These need to eventually be set using configuration file
  tolerance = 1.0e-6;
  max_iteration = 100;

  // Create linear solver
  int t_lsolv = timer->createCategory("Linear Solver");
  timer->start(t_lsolv);
  gridpack::math::LinearSolver isolver(*J);
  isolver.configure(cursor);
  timer->stop(t_lsolv);

  tol = 2.0*tolerance;
  int iter = 0;

  // First iteration
  X->zero(); //might not need to do this
  busIO.header("\nCalling solver\n");
  timer->start(t_lsolv);
  isolver.solve(*PQ, *X);
  timer->stop(t_lsolv);
  tol = X->norm2();
//  busIO.header("\nX values\n");
//  X->print();

  // Exchange new values
  int t_updt = timer->createCategory("Bus Update");
  timer->start(t_updt);
  network->updateBuses();
  timer->stop(t_updt);

  while (real(tol) > tolerance && iter < max_iteration) {
    // Push current values in X vector back into network components
    // Need to implement setValues method in PFBus class in order for this to
    // work
    vMap.mapToBus(X);

    // Exchange data between ghost buses (I don't think we need to exchange data
    // between branches)
    timer->start(t_updt);
    network->updateBuses();
    timer->stop(t_updt);

    // Create new versions of Jacobian and PQ vector
    factory.setMode(RHS);
    vMap.mapToVector(PQ);
//    sprintf(ioBuf,"\nIteration %d Print PQ\n",iter+1);
//    busIO.header(ioBuf);
//    PQ->print();
    factory.setMode(Jacobian);
    jMap.mapToMatrix(J);

    // Create linear solver
    timer->start(t_lsolv);
    gridpack::math::LinearSolver solver(*J);
    solver.configure(cursor);
//    sprintf(ioBuf,"\nIteration %d Print X\n",iter+1);
//    busIO.header(ioBuf);
    X->zero(); //might not need to do this
    solver.solve(*PQ, *X);
    timer->stop(t_lsolv);
//    X->print();

    tol = X->norm2();
    sprintf(ioBuf,"\nIteration %d Tol: %12.6e\n",iter+1,real(tol));
    busIO.header(ioBuf);
    iter++;
  }
#endif

  gridpack::serial_io::SerialBranchIO<PFNetwork> branchIO(128,network);
  branchIO.header("\n   Branch Power Flow\n");
  branchIO.header("\n        Bus 1       Bus 2            P"
                  "                    Q\n");
  branchIO.write();


  busIO.header("\n   Bus Voltages and Phase Angles\n");
  busIO.header("\n   Bus Number      Phase Angle      Voltage Magnitude\n");
  busIO.write();
#endif
  timer->dump();

}
