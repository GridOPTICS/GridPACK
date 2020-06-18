/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_app.cpp
 * @author Shrirang Abhyankar
 * @date   2019-11-21 07:39:01 d3g096
 * 
 * @brief  
 * Example for testing PETSc's implicit solvers
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include "ds_app.hpp"
#include "ds_factory.hpp"
#include "gridpack/include/gridpack.hpp"
#include "gridpack/math/dae_solver.hpp"

namespace gridpack {
namespace dsimplicit {

class DSProblem
  : private gridpack::utility::Uncopyable
{
  // DSFactory
  gridpack::dsimplicit::DSFactory& p_factory;

  //DSNetwork
  boost::shared_ptr<gridpack::dsimplicit::DSNetwork>& p_network;

  gridpack::mapper::BusVectorMap<gridpack::dsimplicit::DSNetwork>& p_VecMapper;
  gridpack::mapper::FullMatrixMap<gridpack::dsimplicit::DSNetwork>& p_MatMapper;

  int p_localsize;

public:
  // Default constructor
  DSProblem(gridpack::dsimplicit::DSFactory& factory,boost::shared_ptr<gridpack::dsimplicit::DSNetwork>& network,gridpack::mapper::BusVectorMap<gridpack::dsimplicit::DSNetwork>& VecMapper, gridpack::mapper::FullMatrixMap<gridpack::dsimplicit::DSNetwork>& MatMapper,int localsize)
  : p_factory(factory), 
    p_network(network), 
    p_VecMapper(VecMapper), 
    p_MatMapper(MatMapper),
    p_localsize(localsize)
  {
    J = p_MatMapper.mapToMatrix();
  }

  // Destructor
  ~DSProblem(void)
  {}

  // A place to build/store the Jacobian
  // We should be able to reuse the Jacobian created for DAESolve.
  // Need to somehow eliminate the need for creating this second Jacobian.
  // The nonlinear solver interface needs to be changed.
  boost::shared_ptr<gridpack::math::Matrix> J;

  // Build the Jacobian for the nonlinear solver at tfaulton or tfaultoff
  void
  operator() (const gridpack::math::Vector& X,gridpack::math::Matrix& J)
  {
    // Push current values in X vector back into network components
    p_factory.setMode(XVECTOBUS);
    p_VecMapper.mapToBus(X);

    // Update ghost buses
    p_network->updateBuses();

    // Evaluate the fault residual Jacobian
    p_factory.setMode(FAULT_EVAL);
    p_MatMapper.mapToMatrix(J);
    
  }

  /// Build the DAE Jacobian
  void operator() (const double& time, 
		   const gridpack::math::Vector& X, 
		   const gridpack::math::Vector& Xdot, 
		   const double& shift, gridpack::math::Matrix& J)
  {
    p_factory.setTSshift(shift);
    // Push current values in X vector back into network components
    p_factory.setMode(XVECTOBUS);
    p_VecMapper.mapToBus(X);

    // Push current values in Xdot vector back into network components
    p_factory.setMode(XDOTVECTOBUS);
    p_VecMapper.mapToBus(Xdot);

    // Update ghost buses
    p_network->updateBuses();

    // Evaluate the DAE Jacobian
    p_factory.setMode(RESIDUAL_EVAL);
    p_MatMapper.mapToMatrix(J);
    //    J.print();
  }

  // Build the residual for the nonlinear solver at tfaulton and tfaultoff
  void
  operator() (const gridpack::math::Vector& X, gridpack::math::Vector& F)
  {
    // Push current values in X vector back into network components
    p_factory.setMode(XVECTOBUS);
    p_VecMapper.mapToBus(X);

    // Update ghost buses
    p_network->updateBuses();

    // Evaluate the residual f(x) - xdot
    p_factory.setMode(FAULT_EVAL);
    p_VecMapper.mapToVector(F);
    F.ready();
    //    F.print();
  }


  /// Build the DAE RHS function
  void operator() (const double& time, 
		   const gridpack::math::Vector& X, const gridpack::math::Vector& Xdot, 
		   gridpack::math::Vector& F)
  {
    // Push current values in X vector back into network components
    p_factory.setMode(XVECTOBUS);
    p_VecMapper.mapToBus(X);

    // Push current values in Xdot vector back into network components
    p_factory.setMode(XDOTVECTOBUS);
    p_VecMapper.mapToBus(Xdot);

    // Update ghost buses
    p_network->updateBuses();

    // Evaluate the residual f(x) - xdot
    p_factory.setMode(RESIDUAL_EVAL);
    p_VecMapper.mapToVector(F);
    F.ready();
    //F.print();
  }
};

// Calling program for implicit dynamic simulation

/**
 * Basic constructor
 */
gridpack::dsimplicit::DSApp::DSApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::dsimplicit::DSApp::~DSApp(void)
{
}

/**
 * Execute application
 * @param argc number of arguments
 * @param argv list of character strings
 */
void gridpack::dsimplicit::DSApp::execute(int argc, char** argv)
{
  int t_total;                                        // Timer ID: Total time
  int t_fileread;                                     // Timer ID: File read
  int t_partition;                                    // Timer ID: Partition
  int t_setup;                                        // Timer ID: Setup
  int t_createMatVec;                                 // Timer ID: Create matrices and vectors
  gridpack::utility::CoarseTimer *timer;              // Coarse timer
  gridpack::parallel::Communicator world;             // World communicator..same as MPI_COMM_WORLD
  gridpack::utility::Configuration *config;           // configuration file reader
  gridpack::utility::Configuration::CursorPtr cursor; // configuration file cursor
  boost::shared_ptr<DSNetwork> network(new DSNetwork(world)); // network object
  gridpack::parser::PTI23_parser<DSNetwork> parser(network);  // parser object for network and dyr file

  std::string netfilename;                            // Name of the file holding network data
  std::string dyrfilename;                            // Name of the file holding dynamic data

  // Create coarse timer instance
  // Timer, as the name suggests, is not a "single" timer but a collection of
  // timers that can be instantiated via the createCategory() method. Each timer
  // category is identified by a unique integer ID.
  timer = gridpack::utility::CoarseTimer::instance();

  // Create timer category for recording the total time taken by the application
  t_total = timer->createCategory("Total Application");
  // Trigger the timer
  timer->start(t_total);

  // Create configuration file reader
  config = gridpack::utility::Configuration::configuration();
  // Set configuration file with the reader
  config->open("input.xml",world);
  // The configuration cursor manages easily reading the lines in the
  // configuration file. The gridpack::utility::Configuration::getCursor()
  // method sets the cursor to the line having the string identifier given
  // to it. The string identifer, usually, is a top level identifier for an
  // application, for example <Dynamic_simulation>. All of its runtime options
  // can be read subsequently using the gridpack::utility::Configuration::get() method.
  cursor = config->getCursor("Configuration.Dynamic_simulation");

  t_fileread = timer->createCategory("Read Files");
  timer->start(t_fileread);
  // Read the name of the network file
  cursor->get("networkConfiguration",&netfilename);
  parser.parse(netfilename.c_str());
  cursor->get("generatorParameters",&dyrfilename);
  parser.parse(dyrfilename.c_str());
  timer->stop(t_fileread);

  /* TO DO: Read Fault events */

  t_partition = timer->createCategory("Partition Network");
  timer->start(t_partition);
  // partition network
  network->partition();
  timer->stop(t_partition);
  
  t_setup = timer->createCategory("Setup");
  // The INIT_X mode initializes the solution vector and the Jacobian matrix. The Jacobian matrix
  // is filled with all zeros.
  timer->start(t_setup);
  // Once the network is partitioned, the data is stored in objects named "DataCollection". Each bus and
  // branch object has a DataCollection object associated with it. The load() method in the factory class
  // loops over all buses and branches in the network and calls the load method in the application bus and branch
  // objects. The bus and branch load method is used to save a copy of the data stored in the DataCollection
  // object to application bus and branch objects for subsequent calculations needed when the solver is
  // called.
  gridpack::dsimplicit::DSFactory factory(network);   // Application factory
  factory.load();

  // set network components using factory
  // the setComponents method sets up the neighbour information for buses and branches. It
  // sets up information such as what buses cover a branch and which branches are incident at
  // a bus
  factory.setComponents();

  
  // Create solution,residual vectors and Jacobian matrix
  t_createMatVec = timer->createCategory("Create Vector and Matrices");
  timer->start(t_createMatVec);
  factory.setMode(INIT_X);
  gridpack::mapper::BusVectorMap<DSNetwork> VecMapper(network); // Mapper for creating matrices
  gridpack::mapper::FullMatrixMap<DSNetwork> MatMapper(network); // Mapper for creating vectors

  boost::shared_ptr<gridpack::math::Vector> X = VecMapper.mapToVector();
#if DEBUG
  X->print();

  boost::shared_ptr<gridpack::math::Vector> R = VecMapper.mapToVector();
  boost::shared_ptr<gridpack::math::Matrix> J = MatMapper.mapToMatrix();
  factory.setMode(RESIDUAL_EVAL);
  MatMapper.mapToMatrix(J);
  R->zero();
  VecMapper.mapToVector(R);
  //  R->print();
  //  J->print();
  exit(1);
#endif
  timer->stop(t_createMatVec);

  // Set up bus data exchange buffers. Need to decide what data needs to be
  // exchanged
  int t_setx = timer->createCategory("Factory: Set Exchange");
  timer->start(t_setx);
  factory.setExchange();
  timer->stop(t_setx);

  // Create bus data exchange
  int t_updt = timer->createCategory("Bus Update");
  timer->start(t_updt);
  network->initBusUpdate();
  timer->stop(t_updt);

  int lsize= X->localSize();
  double maxtime(0.10);
  int maxsteps(10000);
  DSProblem dsprob(factory,network,VecMapper,MatMapper,lsize);

  gridpack::math::DAESolver::JacobianBuilder daejbuilder = boost::ref(dsprob);
  gridpack::math::DAESolver::FunctionBuilder daefbuilder = boost::ref(dsprob);

  gridpack::math::DAESolver daesolver(world, lsize, daejbuilder, daefbuilder);

  // Get simulation time length
  double tmax;
  cursor->get("simulationTime",&tmax);

  // Read fault parameters
  double faultontime(0.1),faultofftime(0.2);
  int    faultbus(9);
  double Gfault(0.0),Bfault(0.0);
  cursor->get("faultontime",&faultontime);
  cursor->get("faultofftime",&faultofftime);
  cursor->get("faultbus",&faultbus);
  cursor->get("Gfault",&Gfault);
  cursor->get("Bfault",&Bfault);

  // Create nonlinear solver for solving the algebraic equations at fault-on/fault-off time instants
  math::NonlinearSolver::JacobianBuilder jbuildf = boost::ref(dsprob);
  math::NonlinearSolver::FunctionBuilder fbuildf = boost::ref(dsprob);

  boost::scoped_ptr<gridpack::math::NonlinearSolver> nlsolver;
  nlsolver.reset(new gridpack::math::NonlinearSolver(*(dsprob.J),jbuildf,fbuildf));
  nlsolver->configure(cursor);
	      
  // Pre-fault time-stepping
  daesolver.configure(cursor);
  daesolver.initialize(0,0.01,*X);
  daesolver.solve(faultontime,maxsteps);

  // Set fault
  printf("Applying a fault on bus %d at t = %3.2f\n",faultbus,faultontime);
  factory.setfault(faultbus,-Gfault,-Bfault);
  // Solve algebraic fault-on equations
  nlsolver->solve(*X);

  // Fault-on time-stepping
  maxsteps = 10000;
  daesolver.initialize(faultontime,0.01,*X);
  daesolver.solve(faultofftime,maxsteps);

  // Remove fault
  printf("Removing fault on bus %d at t = %3.2f\n",faultbus,faultontime);
  factory.setfault(faultbus,Gfault,Bfault);
  // Solve algebraic fault-on equations
  nlsolver->solve(*X);

  // Post-fault time-stepping
  maxsteps = 10000;
  daesolver.initialize(faultofftime,0.01,*X);
  daesolver.solve(tmax,maxsteps);

  timer->stop(t_setup);
  timer->stop(t_total);
  timer->dump();
}

} // namespace dsimplicit
} // namespace gridpack
