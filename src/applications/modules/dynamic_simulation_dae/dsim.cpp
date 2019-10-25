#include <stdio.h>
#include <dsim.hpp>

DSim::DSim(void)
{
  p_isSetUp = 0;
  p_network.reset(new DSimNetwork(p_comm));
}

DSim::DSim(gridpack::parallel::Communicator comm)
{
  p_isSetUp = 0;
  p_comm = comm;
}

void DSim::setconfigurationfile(const char* configfile)
{
  strcpy(p_configfile,configfile);
}

DSim::~DSim(void)
{
  delete(p_factory);
  delete(p_VecMapper);
  delete(p_MatMapper);
  delete(p_daesolver);
  delete(p_nlsolver);
  if(!rank()) printf("DSim: Finished running simulation\n");
}

void DSim::readnetworkdatafromconfig(void)
{
  gridpack::parser::PTI23_parser<DSimNetwork> parser(p_network);
  std::string netfilename; //Network file name
  std::string dyrfilename; // Dynamic data (.dyr) file

  p_config = gridpack::utility::Configuration::configuration();
  p_config->open(p_configfile,p_comm);
  p_configcursor = p_config->getCursor("Configuration.Dynamic_simulation");

  // Start data read timer 
  p_profiler.startdatareadtimer();
  p_configcursor->get("networkConfiguration",&netfilename);
  parser.parse(netfilename.c_str());
  p_configcursor->get("generatorParameters",&dyrfilename);
  parser.parse(dyrfilename.c_str());
  p_profiler.stopdatareadtimer();
  if(!rank()) printf("DSim: Finished Reading data files %s and %s\n",netfilename.c_str(),dyrfilename.c_str());
}

void DSim::setup()
{
  if(p_isSetUp) return;
  p_profiler.startsetuptimer();
  /* Partition network */
  p_network->partition();
  if(!rank()) printf("DSim:Finished partitioning network\n");

  /* Create factory */
  p_factory = new DSimFactory(p_network);
  /* Load data from Data Collection objects to Bus and Branch components */
  p_factory->load();

  /* Set up connectivity information */
  p_factory->setComponents();
  /* Set up ghost/local status */
  p_factory->initialize();
  // Set up bus data exchange buffers.
  p_factory->setExchange();
  if(!rank()) printf("DSim:Finished setting up factory\n");

  // Create bus data exchange
  p_network->initBusUpdate();

  /* Create mappers and vectors, matrices */
  p_VecMapper = new gridpack::mapper::BusVectorMap<DSimNetwork>(p_network);
  p_MatMapper = new gridpack::mapper::FullMatrixMap<DSimNetwork>(p_network);

  //  p_factory->setMode(INIT_X);

  p_X = p_VecMapper->mapToVector();
  p_J = p_MatMapper->mapToMatrix();

  if(!rank()) printf("DSim:Finished setting up mappers\n");

  // Set up solver
  int lsize = p_X->localSize();

  gridpack::math::DAESolver::JacobianBuilder daejacobian = boost::ref(*this);
  gridpack::math::DAESolver::FunctionBuilder daefunction = boost::ref(*this);

  p_daesolver = new gridpack::math::DAESolver(p_comm,lsize,daejacobian,daefunction);

  // Create nonlinear solver for solving the algebraic equations at fault-on/fault-off time instants
  gridpack::math::NonlinearSolver::JacobianBuilder jbuildf = boost::ref(*this);
  gridpack::math::NonlinearSolver::FunctionBuilder fbuildf = boost::ref(*this);

  p_nlsolver = new gridpack::math::NonlinearSolver(*(this->p_J),jbuildf,fbuildf);
  p_nlsolver->configure(p_configcursor);
  p_daesolver->configure(p_configcursor);

  if(!rank()) printf("DSim:Finished setting up DAE solver\n");

  p_isSetUp = 1;
  p_profiler.stopsetuptimer();
  if(!rank()) printf("DSim:Set up completed\n");
}

void DSim::initialize()
{
  p_factory->setMode(INIT_X);
  p_VecMapper->mapToVector(p_X);
  
  p_network->updateBuses();
  //p_X->print();
}

/**
  This routine needs to be updated
  -- Remove the fault stuff and use events instead. Need functionality in
     GridPACK for handling events.
*/
void DSim::solve()
{
  // Get simulation time length
  double tmax;
  int maxsteps(10000);

  p_configcursor->get("simulationTime",&tmax);

  // Read fault parameters
  double faultontime(0.1),faultofftime(0.2);
  int    faultbus(9);
  double Gfault(0.0),Bfault(0.0);
  p_configcursor->get("faultontime",&faultontime);
  p_configcursor->get("faultofftime",&faultofftime);
  p_configcursor->get("faultbus",&faultbus);
  p_configcursor->get("Gfault",&Gfault);
  p_configcursor->get("Bfault",&Bfault);

  // Pre-fault time-stepping
  p_daesolver->initialize(0,0.01,*p_X);
  p_daesolver->solve(faultontime,maxsteps);
  if(!rank()) printf("DSim:Finished pre-fault simulation\n");

  // Set fault
  printf("Applying a fault on bus %d at t = %3.2f\n",faultbus,faultontime);
  p_factory->setfault(faultbus,-Gfault,-Bfault);
  p_factory->setMode(XVECPRETOBUS);
  p_VecMapper->mapToBus(*p_X);
  // Solve algebraic fault-on equations
  p_nlsolver->solve(*p_X);

  // Fault-on time-stepping
  maxsteps = 10000;
  p_daesolver->initialize(faultontime,0.01,*p_X);
  p_daesolver->solve(faultofftime,maxsteps);
  if(!rank()) printf("DSim:Finished fault-on simulation\n");

  // Remove fault
  printf("Removing fault on bus %d at t = %3.2f\n",faultbus,faultontime);
  p_factory->setfault(faultbus,Gfault,Bfault);
  p_factory->setMode(XVECPRETOBUS);
  p_VecMapper->mapToBus(*p_X);
  // Solve algebraic fault-on equations
  p_nlsolver->solve(*p_X);

  // Post-fault time-stepping
  maxsteps = 10000;
  p_daesolver->initialize(faultofftime,0.01,*p_X);
  p_daesolver->solve(tmax,maxsteps);
  if(!rank()) printf("DSim:Finished post-fault simulation\n");
}
