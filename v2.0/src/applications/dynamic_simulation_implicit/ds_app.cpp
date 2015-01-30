/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_app.cpp
 * @author Shrirang Abhyankar
 * @date   2015-01-23 08:25:44 d3g096
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
  double p_maxtime;
  double p_stepsize;
public:
  // Default constructor
  DSProblem(gridpack::dsimplicit::DSFactory& factory,boost::shared_ptr<gridpack::dsimplicit::DSNetwork>& network,gridpack::mapper::BusVectorMap<gridpack::dsimplicit::DSNetwork>& VecMapper, gridpack::mapper::FullMatrixMap<gridpack::dsimplicit::DSNetwork>& MatMapper,
	    int localsize,double maxtime, double stepsize)
  : p_factory(factory), 
    p_network(network), 
    p_VecMapper(VecMapper), 
    p_MatMapper(MatMapper),
    p_localsize(localsize),
    p_maxtime(maxtime),
    p_stepsize(stepsize)
  {}

  // Destructor
  ~DSProblem(void)
  {}

  /// Build a Jacobian
  void operator() (const double& time, 
		   const gridpack::math::Vector& X, 
		   const gridpack::math::Vector& Xdot, 
		   const double& shift, gridpack::math::Matrix& J)
  {
  }

  /// Build the RHS function
  void operator() (const double& time, 
		   const gridpack::math::Vector& X, const gridpack::math::Vector& Xdot, 
		   gridpack::math::Vector& F)
  {
    X.print();
  }

  void solve(const gridpack::parallel::Communicator& comm,
             gridpack::utility::Configuration::CursorPtr conf,double t0,boost::shared_ptr<gridpack::math::Vector> X)
  {
    int maxsteps = 10000;
    gridpack::math::DAEJacobianBuilder jbuilder = boost::ref(*this);
    gridpack::math::DAEFunctionBuilder fbuilder = boost::ref(*this);
    
    gridpack::math::DAESolver solver(comm, p_localsize, jbuilder, fbuilder);
    solver.configure(conf);

    solver.initialize(t0,p_stepsize,*X);

    solver.solve(p_maxtime,maxsteps);
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
  boost::shared_ptr<gridpack::math::Vector> R(X->clone());
  X->print();

  boost::shared_ptr<gridpack::math::Matrix> J = MatMapper.mapToMatrix();
  factory.setMode(RESIDUAL_EVAL);
  MatMapper.mapToMatrix(J);
  R->zero();
  VecMapper.mapToVector(R);
  R->print();
  J->print();
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
  DSProblem dsprob(factory,network,VecMapper,MatMapper,lsize,0.1,0.01);

  dsprob.solve(world,cursor,0,X);

  timer->stop(t_setup);
  timer->stop(t_total);
  timer->dump();
}

} // namespace dsimplicit
} // namespace gridpack
