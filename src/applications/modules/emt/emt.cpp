#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <emt.hpp>

#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/parser/PTI33_parser.hpp"
#include "gridpack/parser/PTI34_parser.hpp"
#include "gridpack/parser/PTI35_parser.hpp"
#include "gridpack/parser/PTI36_parser.hpp"


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
    p_sim->updateprecon = true;
    
    Emt::EventManager::p_handle(nevent, eventidx, t, state);


    p_sim->p_daesolver->restartstep();

    //    p_sim->p_daesolver->reusepreconditioner(-2); // Update preconditioner

  }

  
};

Emt::Emt(void)
  : p_isSetUp(0),
    reuseprecon_nsteps(1),
    updateprecon(false),
    timestep_prev(1.0),
    p_saveoutput(false),
    fp_monitor(NULL),
    emt_network(new EmtNetwork(p_comm))
{}

Emt::Emt(gridpack::parallel::Communicator comm)
  : p_comm(comm),
    p_isSetUp(0),
    p_saveoutput(false),
    timestep_prev(1.0),
    updateprecon(false),
    reuseprecon_nsteps(1),
    fp_monitor(NULL),
    emt_network(new EmtNetwork(p_comm))
{}

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
  if(fp_monitor) fclose(fp_monitor);
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

void Emt::setMonitors(gridpack::utility::Configuration::CursorPtr p_configcursor)
{
  // Set up output
  gridpack::utility::Configuration::CursorPtr list;
  list = p_configcursor->getCursor("Monitors");
  if(list) {
    p_saveoutput = true;
    gridpack::utility::Configuration::ChildElements monitors;
    list->children(monitors);
    int size = monitors.size();
    int idx;
    std::vector<int> gen_buses;
    std::vector<int> buses;
    std::vector<std::string> gen_ids;
    for(idx = 0; idx < size; idx++) {
      if(strcmp(monitors[idx].name.c_str(),"Generator") == 0) {
	int busnum = monitors[idx].cursor->get("bus",0);
	std::vector<int> bus_local_idx;
	std::string raw_gen_id = monitors[idx].cursor->get("id","1");
	gridpack::utility::StringUtils util;
	std::string gen_id = util.clean2Char(raw_gen_id);
	EmtBus *bus;
	BaseEMTGenModel *gen;
	bus_local_idx = emt_network->getLocalBusIndices(busnum);
	// Only the owning rank monitors; ghost copies are never advanced.
	if(bus_local_idx.size() && emt_network->getActiveBus(bus_local_idx[0])) {
	  bus = dynamic_cast<EmtBus*>(emt_network->getBus(bus_local_idx[0]).get());
	  gen = bus->getGenerator(gen_id);
	  if(gen) {
	    monitored_gens.push_back(gen);
	  }
	}						    
      } else if(strcmp(monitors[idx].name.c_str(),"Bus") == 0) {
	int busnum = monitors[idx].cursor->get("bus",0);
	std::vector<int> bus_local_idx;
	EmtBus *bus;
	bus_local_idx = emt_network->getLocalBusIndices(busnum);
	if(bus_local_idx.size() && emt_network->getActiveBus(bus_local_idx[0])) {
	  bus = dynamic_cast<EmtBus*>(emt_network->getBus(bus_local_idx[0]).get());
	  monitored_buses.push_back(bus);
	}
      }
    }
    
    if(p_configcursor->get("MonitorFile",&p_monitorfile)) {
      fp_monitor = fopen(p_monitorfile.c_str(),"w");
    } else {
      if(!rank()) {
	fp_monitor = fopen("emtoutput.csv","w");
      } else {
	std::string outfile = "emtoutput_"+std::to_string(rank())+".csv";
	fp_monitor = fopen(outfile.c_str(),"w");
      }
    }
  }

  if(p_saveoutput) {
    fprintf(fp_monitor,"%s","t ");
    output_string[0] = '\0';
    char *ptr = output_string;
    int len = 0,slen;
    bool write_data;

    for(int i = 0; i < monitored_buses.size(); i++) {
      char buf[128];
      write_data = monitored_buses[i]->serialWrite(buf,128,"header");
      if(write_data) {
	slen = strlen(buf);
	if(len + slen < 2084) snprintf(ptr,slen+1,"%s",buf);
	len += slen;
	ptr += slen;
      }
    }

    for(int i = 0; i < monitored_gens.size(); i++) {
      char buf[128];
      write_data = monitored_gens[i]->serialWrite(buf,128,"header");
      if(write_data) {
	slen = strlen(buf);
	if(len + slen < 2084) snprintf(ptr,slen+1,"%s",buf);
	len += slen;
	ptr += slen;
      }
    }
    fprintf(fp_monitor,"%s",output_string);
    fprintf(fp_monitor,"%s","\n");
  }
}

void Emt::setup()
{  
  if(p_isSetUp) return;
  p_profiler.startsetuptimer();

  p_configcursor = p_config->getCursor("Configuration.Powerflow");

  std::string networkfilename;
    
  if (p_configcursor->get("networkConfiguration_v36",&networkfilename)) {

    gridpack::parser::PTI36_parser<EmtNetwork> parser(emt_network);
    
    p_configcursor = p_config->getCursor("Configuration.EMT");

    // Read generator data
    std::string filename;
    
    filename = p_configcursor->get("generatorParameters","");

    try {
      if(filename.size() > 0) {
        parser.externalParse(filename.c_str());
      } else {
        throw(filename.size());
      }
    }
    catch(int size) {
      std::cout << "Cannot read dynamic data file %s\n" << filename;
    }
  } else if (p_configcursor->get("networkConfiguration_v35",&networkfilename)) {

    gridpack::parser::PTI35_parser<EmtNetwork> parser(emt_network);
    
    p_configcursor = p_config->getCursor("Configuration.EMT");

    // Read generator data
    std::string filename;
    
    filename = p_configcursor->get("generatorParameters","");

    try {
      if(filename.size() > 0) {
	parser.externalParse(filename.c_str());
      } else {
	throw(filename.size());
      }
    }
    catch(int size) {
      std::cout << "Cannot read dynamic data file %s\n" << filename;
    }
  } else if (p_configcursor->get("networkConfiguration_v34",&networkfilename)) {
    gridpack::parser::PTI34_parser<EmtNetwork> parser(emt_network);

        p_configcursor = p_config->getCursor("Configuration.EMT");

    // Read generator data
    std::string filename;
    
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

  } else if (p_configcursor->get("networkConfiguration_v33",&networkfilename)) {
    gridpack::parser::PTI33_parser<EmtNetwork> parser(emt_network);

        p_configcursor = p_config->getCursor("Configuration.EMT");

    // Read generator data
    std::string filename;
    
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

  } else if (p_configcursor->get("networkConfiguration",&networkfilename)) {
    gridpack::parser::PTI23_parser<EmtNetwork> parser(emt_network);

        p_configcursor = p_config->getCursor("Configuration.EMT");

    // Read generator data
    std::string filename;
    
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
  } else {
    printf("No network configuration file specified\n");
    exit(0);
  }

  p_emtmachineintegrationtype = IMPLICIT;
  
  // Get machine integration algorithm
  std::string macinttype;
  macinttype = p_configcursor->get("machineIntegrationType","");

  if(macinttype == "EXPLICIT") p_emtmachineintegrationtype = EXPLICIT;


  /* Create factory */
  p_factory = new EmtFactory(emt_network);

  /* Load data from Data Collection objects to Bus and Branch components */
  p_factory->load();

  /* Set integration type for machines */
  p_factory->setMachineIntegrationType(p_emtmachineintegrationtype);
  
  /* Set up connectivity information */
  p_factory->setComponents();

  // Set up bus data exchange buffers.
  p_factory->setExchange();
  if(!rank()) printf("Emt:Finished setting up factory\n");

  /* Set up buses and branches */
  p_factory->setup();

  /* Read events from configuration file */
  p_factory->readEvents(p_configcursor);

  // Create bus and branch data exchange
  emt_network->initBusUpdate();

  emt_network->initBranchUpdate();

  /* Create mapper for vector */
  p_VecMapper = new gridpack::mapper::GenVectorMap<EmtNetwork,gridpack::RealType,gridpack::math::RealVector>(emt_network);

  /* Set up global locations for buses and branches
     These indices are used during setting values in
     the Jacobian matrix
  */
  p_factory->setGlobalLocations();

  /* Create mapper for matrix */
  p_MatMapper = new gridpack::mapper::GenMatrixMap<EmtNetwork,gridpack::RealType,gridpack::math::RealMatrix>(emt_network);

  /* Jacobian matrix */
  p_J = p_MatMapper->createMatrix();

  if(!rank()) printf("Emt:Finished setting up mappers\n");

  // Set up output monitors
  setMonitors(p_configcursor);


  // Initialize
  initialize(); 
  if(!rank()) printf("Emt: DAE size %d, Jacobian %d x %d\n", p_X->size(), p_J->rows(), p_J->cols());

  // Set up solver
  int lsize = p_X->localSize();

  RealDAESolver::JacobianBuilder daejacobian = boost::ref(*this);
  RealDAESolver::FunctionBuilder daefunction = boost::ref(*this);
  RealDAESolver::PreStepFunction prestepfunction = boost::ref(*this);
  RealDAESolver::PostStepFunction poststepfunction = boost::ref(*this);
  RealDAESolver::EventManagerPtr eman(new EmtEventManager(this));

  p_factory->setEvents(eman,p_VecMapper);

  gridpack::math::RealMatrix* mat_ptr = p_J.get();
  p_daesolver = new RealDAESolver(p_comm,lsize,mat_ptr,daejacobian,daefunction, eman);
  p_daesolver->configure(p_configcursor);
  p_daesolver->preStep(prestepfunction);
  p_daesolver->postStep(poststepfunction);
  

  // Get simulation time length
  double t(0.0),tstep,tmax;

  // Get simulation time frame from input
  p_configcursor->get("timeStep",&tstep);
  p_configcursor->get("simulationTime",&tmax);
  p_simparams.setTimeStep(tstep);
  p_simparams.setFinalTime(tmax);

  p_daesolver->initialize(t,tstep,*(p_X.get()));


  if(!rank()) printf("Emt:Finished setting up DAE solver\n");

  p_initResidualTol = p_configcursor->get("initResidualTol", 1.0e-6);
  p_abortOnInitResidual = p_configcursor->get("abortOnInitResidual", false);
  if(p_configcursor->get("checkInitialConditions", true)) checkInitialResidual();

  p_isSetUp = 1;
  p_profiler.stopsetuptimer();
  if(!rank()) printf("Emt:Set up completed\n");
}

// Consistent-initialization check. Rows whose residual changes between
// Xdot=0 and Xdot=1 are differential; those with a large required derivative
// (|F0/mass| >= 1) are sinusoidal abc states. Each aligned abc triple gets its
// exact derivative from the phase values (x=A sin(wt+phi), xdot=wA cos(wt+phi));
// all other derivatives are zero. F(0, X0, Xdot) must then vanish on every row
// if X0 is a sinusoidal steady state of the EMT network.
void Emt::checkInitialResidual()
{
  // Zero-length pre-step so explicitly integrated models have computed the
  // internals (field voltage/current, block outputs) the residual reads.
  p_factory->preStep(0.0, 0.0);
  boost::scoped_ptr<gridpack::math::RealVector> xdot(p_X->clone());
  boost::scoped_ptr<gridpack::math::RealVector> f(p_X->clone());
  int lo, hi;
  p_X->localIndexRange(lo, hi);
  int n = hi - lo, i, j, off, nv;
  std::vector<double> x0(n + 1), v0(n + 1), v1(n + 1), v2(n + 1), xd(n + 1, 0.0);
  p_X->getElementRange(lo, hi, &x0[0]);
  xdot->zero(); xdot->ready();
  (*this)(0.0, *p_X, *xdot, *f);
  f->getElementRange(lo, hi, &v0[0]);
  xdot->fill(1.0); xdot->ready();
  (*this)(0.0, *p_X, *xdot, *f);
  f->getElementRange(lo, hi, &v1[0]);

  // Local row -> owning bus/branch and variable index.
  std::vector<int> isBus(n, -1), ownerIdx(n, -1), ownerVar(n, -1);
  int nb = emt_network->numBuses(), nbr = emt_network->numBranches();
  for(i = 0; i < nb; i++) {
    p_VecMapper->getLocalBusOffset(i, &off, &nv);
    for(j = 0; j < nv; j++) {
      if(off + j < 0 || off + j >= n) continue;
      isBus[off+j] = 1; ownerIdx[off+j] = i; ownerVar[off+j] = j;
    }
  }
  for(i = 0; i < nbr; i++) {
    p_VecMapper->getLocalBranchOffset(i, &off, &nv);
    for(j = 0; j < nv; j++) {
      if(off + j < 0 || off + j >= n) continue;
      isBus[off+j] = 0; ownerIdx[off+j] = i; ownerVar[off+j] = j;
    }
  }

  // Sinusoidal derivatives for aligned abc triples within one owner.
  // A phase can have zero derivative at t=0, so a triple counts as sinusoidal
  // when all three rows are differential and any one has |xdot| >= 1.
  std::vector<char> sinus(n, 0), diff(n, 0);
  for(i = 0; i < n; i++) {
    double m = v0[i] - v1[i];
    if(fabs(m) > 1e-12) { diff[i] = 1; if(fabs(v0[i] / m) >= 1.0) sinus[i] = 1; }
  }
  const double s3 = sqrt(3.0);
  int ntriples = 0;
  for(i = 0; i + 2 < n; i++) {
    if(!(diff[i] && diff[i+1] && diff[i+2])) continue;
    if(!(sinus[i] || sinus[i+1] || sinus[i+2])) continue;
    if(isBus[i] != isBus[i+2] || ownerIdx[i] != ownerIdx[i+2] ||
       ownerVar[i+2] != ownerVar[i] + 2) continue;
    sinus[i] = sinus[i+1] = sinus[i+2] = 1;
    ntriples++;
    double asin_ = x0[i], acos_ = (x0[i+2] - x0[i+1]) / s3;
    xd[i]   = OMEGA_S * acos_;
    xd[i+1] = OMEGA_S * (-0.5 * acos_ + 0.5 * s3 * asin_);
    xd[i+2] = OMEGA_S * (-0.5 * acos_ - 0.5 * s3 * asin_);
    i += 2;
  }
  if(n > 0) xdot->setElementRange(lo, hi, &xd[0]);
  xdot->ready();
  (*this)(0.0, *p_X, *xdot, *f);
  f->getElementRange(lo, hi, &v2[0]);
  // Leave the network holding X0 with zero derivatives.
  xdot->zero(); xdot->ready();
  (*this)(0.0, *p_X, *xdot, *f);

  // Optional full row map (<dumpVariableMap>true</dumpVariableMap>), for
  // relating PETSc row numbers to network variables.
  if(p_configcursor->get("dumpVariableMap", false)) {
    for(i = 0; i < n; i++) {
      std::string label = "unmapped";
      if(isBus[i] == 1)
        label = dynamic_cast<EmtBus*>(emt_network->getBus(ownerIdx[i]).get())->describeVariable(ownerVar[i]);
      else if(isBus[i] == 0) {
        EmtBranch *br = dynamic_cast<EmtBranch*>(emt_network->getBranch(ownerIdx[i]).get());
        char b[96];
        snprintf(b, sizeof(b), "branch %d-%d x[%d]", br->getBus1OriginalIndex(), br->getBus2OriginalIndex(), ownerVar[i]);
        label = b;
      }
      printf("  row %6d  x=%13.6e  %s%s\n", lo + i, x0[i], label.c_str(),
             sinus[i] ? " [abc]" : (diff[i] ? " [state]" : " [algebraic]"));
    }
  }
  double maxres = 0.0;
  int nbad = 0;
  std::vector<std::string> lines;
  for(i = 0; i < n; i++) {
    double val = fabs(v2[i]);
    if(val > maxres) maxres = val;
    if(val <= p_initResidualTol) continue;
    nbad++;
    std::string label = "unmapped row";
    if(isBus[i] == 1) {
      label = dynamic_cast<EmtBus*>(emt_network->getBus(ownerIdx[i]).get())
                ->describeVariable(ownerVar[i]);
    } else if(isBus[i] == 0) {
      EmtBranch *br = dynamic_cast<EmtBranch*>(emt_network->getBranch(ownerIdx[i]).get());
      char b[96];
      snprintf(b, sizeof(b), "branch %d-%d x[%d]",
               br->getBus1OriginalIndex(), br->getBus2OriginalIndex(), ownerVar[i]);
      label = b;
    }
    char buf[192];
    snprintf(buf, sizeof(buf), "  residual %12.4e  (%s)  %s", v2[i],
             sinus[i] ? "abc" : (fabs(v0[i]-v1[i]) > 1e-12 ? "state" : "algebraic"),
             label.c_str());
    lines.push_back(buf);
  }
  p_comm.max(&maxres, 1);
  int cnt[2] = { nbad, ntriples };
  p_comm.sum(cnt, 2);
  if(!rank()) {
    printf("Emt: initial-condition check: max residual %.3e with sinusoidal derivatives "
           "(%d abc triples), %d rows above tol %.1e\n", maxres, cnt[1], cnt[0], p_initResidualTol);
  }
  for(int r = 0; r < size(); r++) {   // rank by rank to avoid interleaving
    if(r == rank()) {
      for(size_t k = 0; k < lines.size() && k < 40; k++) {
        if(size() > 1) printf("  [rank %d]", r);
        printf("%s\n", lines[k].c_str());
      }
      if(lines.size() > 40) printf("  ... %zu more rows\n", lines.size() - 40);
      fflush(stdout);
    }
    p_comm.barrier();
  }
  if(cnt[0] > 0) {
    if(p_abortOnInitResidual) {
      if(!rank()) printf("Emt: ERROR inconsistent initial conditions; aborting (abortOnInitResidual=true)\n");
      MPI_Abort(static_cast<MPI_Comm>(p_comm), 1);
    } else if(!rank()) {
      printf("Emt: WARNING initial conditions are not an equilibrium; the case will drift "
             "(set <abortOnInitResidual>true</abortOnInitResidual> to stop here)\n");
    }
  }
}

void Emt::initialize()
{
  p_factory->setMode(INIT_X);
  p_factory->setTime(0.0);
  p_X = p_VecMapper->mapToVector();
  p_X->ready();
  
  emt_network->updateBuses();
  emt_network->updateBranches();

  save_output();
  //p_X->print();
}

void Emt::solve()
{
  int maxsteps(1000000);
  double final_time;
  p_profiler.startsolvetimer();
  p_simparams.getFinalTime(&final_time);
  p_daesolver->solve(final_time,maxsteps);
  p_profiler.stopsolvetimer();
}

void Emt::save_output(const double& time)
{
  if(p_saveoutput) {
    fprintf(fp_monitor,"%6.5f",time);
    output_string[0] = '\0';
    char *ptr = output_string;
    int len = 0,slen;
    bool write_data;
    
    for(int i = 0; i < monitored_buses.size(); i++) {
      char buf[128];
      write_data = monitored_buses[i]->serialWrite(buf,128,"monitor");
      if(write_data) {
	slen = strlen(buf);
	if(len + slen < 2084) snprintf(ptr,slen+1,"%s",buf);
	len += slen;
	ptr += slen;
      }
    }
    
    for(int i = 0; i < monitored_gens.size(); i++) {
      char buf[128];
      write_data = monitored_gens[i]->serialWrite(buf,128,"monitor");
      if(write_data) {
	slen = strlen(buf);
	if(len + slen < 2084) snprintf(ptr,slen+1,"%s",buf);
	len += slen;
	ptr += slen;
      }
    }
    fprintf(fp_monitor,"%s",output_string);
    fprintf(fp_monitor,"%s","\n");
  }
}

