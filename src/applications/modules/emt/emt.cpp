#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <emt.hpp>

// -------------------------------------------------------------
//  class EmtEventManager
// -------------------------------------------------------------
class EmtEventManager
  : public Emt::EventManager
{
public:

  /// Default constructor.
  EmtEventManager(Emt *sim)
    : Emt::EventManager(), p_sim(sim)
  {}

  /// Destructor
  ~EmtEventManager(void)
  {}

protected:

  /// The simulation
  Emt *p_sim;

  /// Handle triggered events (specialized)
  virtual void p_handle(const int nevent,
                        int *eventidx, const double t,
                        Emt::VectorType& state)
  {

    //    p_sim->p_X->equate(state);
    //    p_sim->p_factory->setMode(XVECTOBUS);
    //    p_sim->p_VecMapper->mapToBus(*(p_sim->p_X));
  
    Emt::EventManager::p_handle(nevent, eventidx, t, state);

  }

  
};


// -------------------------------------------------------------
// class EmtTimedFaultEvent
// 
// This manages 2 solver events: when the fault goes on and when it
// goes off.
// -------------------------------------------------------------
class EmtTimedFaultEvent
  : public Emt::Event
{
public:

  /// Default constructor.
  EmtTimedFaultEvent(Emt *sim,
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
  ~EmtTimedFaultEvent(void)
  {}

protected:

  Emt *p_sim;
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
  } else if (triggered[1]) {
    // Remove fault
  }
}

};


Emt::Emt(void)
{
  p_isSetUp = 0;
  emt_network.reset(new EmtNetwork(p_comm));
}

Emt::Emt(gridpack::parallel::Communicator comm)
{
  p_isSetUp = 0;
  p_comm = comm;
}

void Emt::setconfigurationfile(const char* configfile)
{
  // FIXME: Really should use strncpy(), but what is wrong with using
  // std::string?
  strcpy(p_configfile,configfile);
}

Emt::~Emt(void)
{
  delete(p_factory);
  delete(p_VecMapper);
  delete(p_MatMapper);
  delete(p_pfapp);
  delete(p_daesolver);
  if(!rank()) printf("Emt: Finished running simulation\n");
}

/**
 * Transfer data from power flow to EMT simulattion
 * @param pf_network power flow network
 * @param emt_network EMT simulation network
 */
void Emt::transferPFtoEMT(
    boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network,
    boost::shared_ptr<EmtNetwork>
    emt_network)
{
  int numBus = pf_network->numBuses();
  int i;
  gridpack::component::DataCollection *pfData;
  gridpack::component::DataCollection *emtData;
  double rval;
  for (i=0; i<numBus; i++) {
    pfData = pf_network->getBusData(i).get();
    emtData = emt_network->getBusData(i).get();
    pfData->getValue("BUS_PF_VMAG",&rval);
    emtData->setValue(BUS_VOLTAGE_MAG,rval);
    ///printf("Step0 bus%d mag = %f\n", i+1, rval);
    pfData->getValue("BUS_PF_VANG",&rval);
    emtData->setValue(BUS_VOLTAGE_ANG,rval);
    int ngen = 0;
    if (pfData->getValue(GENERATOR_NUMBER, &ngen)) {
      int j;
      for (j=0; j<ngen; j++) {
        pfData->getValue("GENERATOR_PF_PGEN",&rval,j);
        emtData->setValue(GENERATOR_PG,rval,j);
        //printf("save PGEN: %f\n", rval);
        pfData->getValue("GENERATOR_PF_QGEN",&rval,j);
        emtData->setValue(GENERATOR_QG,rval,j);
        //printf("save QGEN: %f\n", rval);
      }
    }
  }
}


void Emt::solvepowerflow(void)
{
  p_config = gridpack::utility::Configuration::configuration();
  p_config->open(p_configfile,p_comm);

  p_configcursor = p_config->getCursor("Configuration.Powerflow");
  bool useNonLinear = false;
  useNonLinear = p_configcursor->get("UseNonLinear", useNonLinear);

  boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network(new gridpack::powerflow::PFNetwork(p_comm));
  
  p_pfapp = new gridpack::powerflow::PFAppModule;
  p_pfapp->readNetwork(pf_network, p_config);
  p_pfapp->initialize();
  if (useNonLinear) {
    p_pfapp->nl_solve();
  } else {
    p_pfapp->solve();
  }
  p_pfapp->write();
  p_pfapp->saveData();

  pf_network->clone<EmtBus,EmtBranch>(emt_network);

  transferPFtoEMT(pf_network,emt_network);
}

void Emt::readnetworkdatafromconfig(void)
{
  gridpack::parser::PTI23_parser<EmtNetwork> parser(emt_network);
  std::string netfilename; //Network file name
  std::string dyrfilename; // Dynamic data (.dyr) file

  p_config = gridpack::utility::Configuration::configuration();
  p_config->open(p_configfile,p_comm);

  // Start data read timer 
  p_profiler.startdatareadtimer();
  p_configcursor->get("networkConfiguration",&netfilename);
  parser.parse(netfilename.c_str());
  p_configcursor->get("generatorParameters",&dyrfilename);
  parser.parse(dyrfilename.c_str());
  p_profiler.stopdatareadtimer();
  if(!rank()) printf("Emt: Finished Reading data files %s and %s\n",netfilename.c_str(),dyrfilename.c_str());
}

void Emt::setup()
{
  if(p_isSetUp) return;
  p_profiler.startsetuptimer();

  p_configcursor = p_config->getCursor("Configuration.Dynamic_simulation");

  // Read generator data
  std::string filename;
  gridpack::parser::PTI23_parser<EmtNetwork> parser(emt_network);

  filename = p_configcursor->get("generatorParameters","");
  try {
    if(filename.size() > 0) parser.externalParse(filename.c_str());
    else {
      throw(filename.size());
    }
  }
  catch(int size) {
    std::cout << "Cannot read dynamic data file %s\n" << filename;
  }
    
  /* Create factory */
  p_factory = new EmtFactory(emt_network);

  /* Load data from Data Collection objects to Bus and Branch components */
  p_factory->load();

  /* Set up connectivity information */
  p_factory->setComponents();
  /* Set up ghost/local status */
  p_factory->initialize();
  // Set up bus data exchange buffers.
  p_factory->setExchange();
  if(!rank()) printf("Emt:Finished setting up factory\n");

  // Create bus data exchange
  emt_network->initBusUpdate();

  /* Create mappers and vectors, matrices */
  p_factory->setMode(INIT_X);

  p_VecMapper = new gridpack::mapper::GenVectorMap<EmtNetwork>(emt_network);

  p_X = p_VecMapper->mapToVector();

  p_MatMapper = new gridpack::mapper::GenMatrixMap<EmtNetwork>(emt_network);

  p_J = p_MatMapper->mapToMatrix();

  if(!rank()) printf("Emt:Finished setting up mappers\n");

  // Initialize
  initialize(); 

  // Set up solver
  int lsize = p_X->localSize();

  DAESolver::JacobianBuilder daejacobian = boost::ref(*this);
  DAESolver::FunctionBuilder daefunction = boost::ref(*this);
  DAESolver::EventManagerPtr eman(new EmtEventManager(this));

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

  EventPtr e(new EmtTimedFaultEvent(this, faultontime, faultofftime,
                                     faultbus, Gfault, Bfault));

  // Event manager must be populated before DAE solver construction
  eman->add(e);

  gridpack::math::Matrix* mat_ptr = p_J.get();
  p_daesolver = new DAESolver(p_comm,lsize,mat_ptr,daejacobian,daefunction, eman);
  p_daesolver->configure(p_configcursor);

  // Get simulation time length
  double t(0.0),tstep,tmax;

  // Get simulation time frame from input
  p_configcursor->get("timeStep",&tstep);
  p_configcursor->get("simulationTime",&tmax);
  p_simparams.setTimeStep(tstep);
  p_simparams.setFinalTime(tmax);

  p_daesolver->initialize(t,tstep,*(p_X.get()));


  if(!rank()) printf("Emt:Finished setting up DAE solver\n");

  p_isSetUp = 1;
  p_profiler.stopsetuptimer();
  if(!rank()) printf("Emt:Set up completed\n");
}

void Emt::initialize()
{
  //  p_factory->setMode(INIT_X);
  p_VecMapper->mapToVector(p_X);
  p_X->ready();
  
  emt_network->updateBuses();
  //p_X->print();
}

/**
  This routine needs to be updated
  -- Remove the fault stuff and use events instead. Need functionality in
     GridPACK for handling events.
*/
void Emt::solve()
{
  int maxsteps(10000);
  double final_time;
  p_profiler.startsolvetimer();
  p_simparams.getFinalTime(&final_time);
  p_daesolver->solve(final_time,maxsteps);
  p_profiler.stopsolvetimer();
}

