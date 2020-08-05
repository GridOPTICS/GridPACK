

#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <dsim.hpp>

// -------------------------------------------------------------
//  class DSimEventManager
// -------------------------------------------------------------
class DSimEventManager
  : public DSim::EventManager
{
public:

  /// Default constructor.
  DSimEventManager(DSim *sim)
    : DSim::EventManager(), p_sim(sim)
  {}

  /// Destructor
  ~DSimEventManager(void)
  {}

protected:

  /// The simulation
  DSim *p_sim;

  /// Handle triggered events (specialized)
  virtual void p_handle(const int nevent,
                        int *eventidx, const double t,
                        DSim::VectorType& state)
  {
    p_sim->p_resolve = false;

    //    p_sim->p_X->equate(state);
    //    p_sim->p_factory->setMode(XVECTOBUS);
    //    p_sim->p_VecMapper->mapToBus(*(p_sim->p_X));
  
    DSim::EventManager::p_handle(nevent, eventidx, t, state);

    // FIXME: All reduce on p_resolve needed here
    if (p_sim->p_resolve) {
      // get the current DAE solution onto the network
      p_sim->p_X->equate(state);
      p_sim->p_factory->setMode(XVECPRETOBUS);
      p_sim->p_VecMapper->mapToBus(*(p_sim->p_X));
      // Solve algebraic equations 
      p_sim->p_nlsolver->solve(*(p_sim->p_X));

      // Push the updated solution back to network
      p_sim->p_factory->setMode(XVECTOBUS);
      p_sim->p_VecMapper->mapToBus(*(p_sim->p_X));
      p_sim->p_X->ready();
      p_sim->p_factory->resetEventFlags();
    }
  }

  
};


// -------------------------------------------------------------
// class DSimTimedFaultEvent
// 
// This manages 2 solver events: when the fault goes on and when it
// goes off.
// -------------------------------------------------------------
class DSimTimedFaultEvent
  : public DSim::Event
{
public:

  /// Default constructor.
  DSimTimedFaultEvent(DSim *sim,
                      const double& ton,
                      const double& toff,
                      const int& busno,
                      const double& Gfault,
                      const double& Bfault)
    : gridpack::math::DAESolver::Event(2),
      p_sim(sim), p_ton(ton), p_toff(toff), p_bus(busno),
      p_Gfault(Gfault), p_Bfault(Bfault)
  {
    // A fault requires that the DAE solver be reset. 
    std::fill(p_term.begin(), p_term.end(), false);
    
    // The event occurs when the values go from positive to negative.
    std::fill(p_dir.begin(), p_dir.end(),
              gridpack::math::CrossZeroNegative);
  }

  /// Destructor
  ~DSimTimedFaultEvent(void)
  {}

protected:

  DSim *p_sim;
  double p_ton, p_toff;
  int p_bus;
  double p_Gfault;
  double p_Bfault;

  void p_update(const double& t, gridpack::ComplexType *state)
  {
    p_current[0] = p_ton - t;
    p_current[1] = p_toff - t;
  }

  void p_handle(const bool *triggered, const double& t,
                gridpack::ComplexType *state)
  {
  // FIXME: Is setfault local (to the process w/ p_bus or collective?
  // The Event class is *local*, so if collective operations are
  // needed they should be handled by an EventManager
    
  if (triggered[0]) {
    // Set fault
    if(!p_sim->rank())printf("Applying a fault on bus %d at t = %3.2f\n", p_bus, t);
    p_sim->p_factory->setfault(p_bus, p_Gfault, -p_Bfault);
    p_sim->p_resolve = true;
  } else if (triggered[1]) {
    // Remove fault
    if(!p_sim->rank()) printf("Removing fault on bus %d at t = %3.2f\n", p_bus, t);
    p_sim->p_factory->setfault(p_bus,-p_Gfault,p_Bfault);
    p_sim->p_resolve = true;
  }
}


};


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
  // FIXME: Really should use strncpy(), but what is wrong with using
  // std::string?
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

  p_X = p_VecMapper->mapToVector();
  p_J = p_MatMapper->mapToMatrix();

  if(!rank()) printf("DSim:Finished setting up mappers\n");

  // Initialize
  initialize(); 

  // Set up solver
  int lsize = p_X->localSize();

  DAESolver::JacobianBuilder daejacobian = boost::ref(*this);
  DAESolver::FunctionBuilder daefunction = boost::ref(*this);
  DAESolver::EventManagerPtr eman(new DSimEventManager(this));

  p_factory->setEvents(eman,p_VecMapper);
  // Read fault parameters, Set up fault events

  double faultontime(0.1),faultofftime(0.2);
  int    faultbus(9);
  double Gfault(0.0),Bfault(0.0);
  p_configcursor->get("faultontime",&faultontime);
  p_configcursor->get("faultofftime",&faultofftime);
  p_configcursor->get("faultbus",&faultbus);
  p_configcursor->get("Gfault",&Gfault);
  p_configcursor->get("Bfault",&Bfault);

  EventPtr e(new DSimTimedFaultEvent(this, faultontime, faultofftime,
                                     faultbus, Gfault, Bfault));

  // Event manager must be populated before DAE solver construction
  eman->add(e);

  gridpack::math::Matrix* mat_ptr = p_J.get();
  p_daesolver = new DAESolver(p_comm,lsize,mat_ptr,daejacobian,daefunction, eman);
  p_daesolver->configure(p_configcursor);

  // Create nonlinear solver for solving the algebraic equations at fault-on/fault-off time instants
  gridpack::math::NonlinearSolver::JacobianBuilder jbuildf = boost::ref(*this);
  gridpack::math::NonlinearSolver::FunctionBuilder fbuildf = boost::ref(*this);

  p_nlsolver = new gridpack::math::NonlinearSolver(*(this->p_J),jbuildf,fbuildf);
  p_nlsolver->configure(p_configcursor);

  if(!rank()) printf("DSim:Finished setting up DAE solver\n");

  p_isSetUp = 1;
  p_profiler.stopsetuptimer();
  if(!rank()) printf("DSim:Set up completed\n");
}

void DSim::initialize()
{
  p_factory->setMode(INIT_X);
  p_VecMapper->mapToVector(p_X);
  p_X->ready();
  
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
  double t(0.0), tstep, tmax;
  int maxsteps(10000);

  // Get simulation time frame from input
  p_configcursor->get("timeStep",&tstep);
  p_configcursor->get("simulationTime",&tmax);

  p_daesolver->initialize(t,tstep,*(p_X.get()));
  p_daesolver->solve(tmax,maxsteps);
  /*
  while (t < tmax) {
    double tout(tmax);
    int maxsteps(10000);
   
    p_daesolver->initialize(t, tstep, *p_X);
    p_daesolver->solve(tout, maxsteps);
    std::cout << "Time requested = " << tmax << ", "
              << "actual time = " << tout << ", "
              << "Steps = " << maxsteps << std::endl;
    t = tout;
  } 
  */ 
}

