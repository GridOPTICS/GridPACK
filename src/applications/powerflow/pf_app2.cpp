// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   pf_app2.cpp
 * @author William A. Perkins
 * @date   2014-01-28 11:30:04 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <ga++.h>

#include "pf_app2.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/math/newton_raphson_solver.hpp"
#include "gridpack/math/nonlinear_solver.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/utilities/uncopyable.hpp"


namespace gridpack {
namespace powerflow {

/// A helper functor for the powerflow solver.
/**
 * This is a utility functor that provides functions that build the
 * Jacobian and RHS from the network.  
 * 
 */
struct PFSolverHelper 
  : private utility::Uncopyable
{

  /// The powerflow factory 
  PFFactory& p_factory;

  /// The powerflow network controlled by ::p_factory
  boost::shared_ptr<PFNetwork> p_network;

  /// A place to build/store the Jacobian
  boost::shared_ptr<math::Matrix> J;

  /// The network state estimate from previous solver iteration
  /**
   * See ::update() for why this is necessary.
   * 
   */
  boost::shared_ptr<math::Vector> Xold;

  /// The current network state estimate
  /**
   * This vector provides a space for the nonlinear solve to store the
   * current solution estimate.  It should be filled with the initial
   * condition, handed to the nonlinear solver, and not changed
   * afterward.
   * 
   */
  boost::shared_ptr<math::Vector> X;

  /// The difference between the current and previous estimate.
  /**
   * See ::update() for why this is necessary.
   * 
   */
  boost::shared_ptr<math::Vector> Xdelta;

  /// Constructor
  /** 
   * The current network state is gathered from the network.
   * 
   * @param factory powerflow factory 
   * @param network network controlled by @c factory
   * 
   * @return 
   */
  PFSolverHelper(PFFactory& factory, boost::shared_ptr<PFNetwork> network)
    : p_factory(factory), p_network(network), Xold(), Xdelta()
  {
    p_factory.setMode(State);
    mapper::BusVectorMap<PFNetwork> vMap(p_network);
    Xold = vMap.mapToVector();
    // Xold->print();
    X.reset(Xold->clone());
    Xdelta.reset(Xold->clone());
    Xdelta->zero();
    p_factory.setMode(Jacobian);
    mapper::FullMatrixMap<PFNetwork> jMap(p_network);
    J = jMap.mapToMatrix();
  }
  
  /// Push the current estimated state back onto the network 
  /** 
   * The network state (voltage, phase) is updated with the current
   * estimate from the nonlinear solver.
   * 
   * FIXME: The problem here is that ...::mapToBus expects the
   * @e CHANGE (old - current)? in state variables. IMHO, there should
   * be a way to set the state variables directly.  The solver should
   * be responsible for making the @e entire estimate.
   *
   * @param Xcur current state estimate from the solver
   */
  void
  update(const math::Vector& Xcur)
  {
    
    Xdelta->equate(Xcur);
    Xdelta->scale(-1.0);
    Xdelta->add(*Xold);
    double snorm(abs(Xdelta->norm2()));
    if (Xdelta->processor_rank() == 0) {
      std::cout << "PFSolverHelper::update(): solution residual: " << snorm << std::endl;
    }

    // Xdelta->print();
    p_factory.setMode(RHS);
    mapper::BusVectorMap<PFNetwork> vMap(p_network);
    vMap.mapToBus(Xdelta);
    Xold->equate(Xcur);
    
    // Exchange data between ghost buses (I don't think we need to exchange data
    // between branches)
    p_network->updateBuses();
    
  }
  
  /// Build the Jacobian Matrix
  /** 
   * This is called by the nonlinear solver each iteration to build
   * the Jacobian from the current network state.  
   *
   * @param Xcur current state estimate
   * @param J Jacobian 
   */
  void
  operator() (const math::Vector& Xcur, math::Matrix& theJ)
  {
    // In both the Netwon-Raphson and PETSc nonlinear solver (some
    // methods) implementations, the RHS function builder is called
    // before this, so we may be able to count on the current solution
    // being on the netork when here.

    // X.print();
    // update(Xcur);
    
    // Set to build Jacobian
    p_factory.setMode(Jacobian);
    mapper::FullMatrixMap<PFNetwork> jMap(p_network);
    
    // build the Jacobian
    jMap.mapToMatrix(theJ);
  }
  
  /// Build the RHS function vector
  /** 
   * This is called by the nonlinear solver each iteration to build
   * the RHS vector from the current network state.  This is also
   * responsible for updating the network state with the current
   * solution estimate, which assumes that it is called only once per
   * solver iteration.  This may not be true if certain methods are
   * used or if a finite difference Jacobian is computed by the
   * solver.
   * 
   * @param Xcur current state estimate
   * @param PQ computed RHS vector
   */
  void
  operator() (const math::Vector& Xcur, math::Vector& PQ)
  {
    // In both the Netwon-Raphson and PETSc nonlinear solver
    // implementations, this is called before the Jacobian builder, so
    // we may only need to map the solution back on to the network here.
    
    // X.print();
    update(Xcur);
    
    // set to build RHS vector
    p_factory.setMode(RHS);
    mapper::BusVectorMap<PFNetwork> vMap(p_network);
    
    // build the RHS vector
    vMap.mapToVector(PQ);
  }
};

// -------------------------------------------------------------
//  class PFApp2
// -------------------------------------------------------------

// -------------------------------------------------------------
// PFApp2:: constructors / destructor
// -------------------------------------------------------------
PFApp2::PFApp2()
{
  
}

PFApp2::~PFApp2(void)
{
}

// -------------------------------------------------------------
// PFApp2::execute
// -------------------------------------------------------------
void
PFApp2::execute(void)
{
  parallel::Communicator world;
  boost::shared_ptr<PFNetwork> network(new PFNetwork(world));

  // read configuration file
  utility::Configuration *config = 
    utility::Configuration::configuration();
  config->open("input.xml", world);
  utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Powerflow");
  std::string filename = cursor->get("networkConfiguration",
      "No network configuration specified");

  // load input file
  utility::CoarseTimer *timer =
    utility::CoarseTimer::instance();
  int t_pti = timer->createCategory("PTI Parser");
  timer->start(t_pti);
  parser::PTI23_parser<PFNetwork> parser(network);
  parser.parse(filename.c_str());
  timer->stop(t_pti);

  // partition network
  int t_part = timer->createCategory("Partition");
  timer->start(t_part);
  network->partition();
  timer->stop(t_part);

  // create factory
  PFFactory factory(network);
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

  // make Sbus components to create S vector
  factory.setSBus();

  // Solve the problem

  bool useNewton(false);
  useNewton = cursor->get("UseNewton", useNewton);

  int t_lsolv = timer->createCategory("Solver");

  timer->start(t_lsolv);

  PFSolverHelper helper(factory, network);
  math::JacobianBuilder jbuildf = boost::ref(helper);
  math::FunctionBuilder fbuildf = boost::ref(helper);

  boost::scoped_ptr<math::NonlinearSolverInterface> solver;
  if (useNewton) {
    math::NewtonRaphsonSolver *tmpsolver = 
      new math::NewtonRaphsonSolver(*(helper.J), jbuildf, fbuildf);
    tmpsolver->tolerance(1.0e-08);
    tmpsolver->maximumIterations(50);
    solver.reset(tmpsolver);
  } else {
    solver.reset(new math::NonlinearSolver(*(helper.J), jbuildf, fbuildf));
  }

  try {
    solver->configure(cursor);
    solver->solve(*helper.X);
    helper.update(*helper.X);
  } catch (const Exception& e) {
    std::cerr << e.what() << std::endl;
  }

  timer->stop(t_lsolv);

  // Report the results

  serial_io::SerialBusIO<PFNetwork> busIO(128,network);
  serial_io::SerialBranchIO<PFNetwork> branchIO(128,network);

  branchIO.header("\n   Branch Power Flow\n");
  branchIO.header("\n        Bus 1       Bus 2            P"
                  "                    Q\n");
  branchIO.write();


  busIO.header("\n   Bus Voltages and Phase Angles\n");
  busIO.header("\n   Bus Number      Phase Angle      Voltage Magnitude\n");
  busIO.write();

  timer->dump();
}


} // namespace powerflow
} // namespace gridpack


// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  gridpack::parallel::Communicator world;
  gridpack::math::Initialize();
  GA_Initialize();
  MA_init(MT_C_CHAR, 1024*1024, 1024*1024);

  gridpack::powerflow::PFApp2 app;
  app.execute();

  GA_Terminate();

  gridpack::math::Finalize();
}

