/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wind_driver.cpp
 * @author Bruce Palmer
 * @date   February 10, 2014
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#include <algorithm>
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/dynamic_simulation_full_y/dsf_app_module.hpp"
#include "wind_driver.hpp"

// Sets up multiple communicators so that individual contingency calculations
// can be run concurrently

/**
 * Basic constructor
 * @param comm communicator used for analysis
 * @param nwatch number of generator variables being watched
 * @param nconf number of scenarios
 * @param nsteps number of timesteps being stored
 */
gridpack::contingency_analysis::QuantileAnalysis::QuantileAnalysis(
    gridpack::parallel::Communicator comm, int nwatch, int nconf, int nsteps)
{
  int dims[3];
  int three = 3;
  int chunk[3];
  p_nwatch = nwatch;
  p_nconf = nconf;
  p_nsteps = nsteps;
  dims[0] = nconf;
  dims[1] = nwatch;
  dims[2] = nsteps;
  // Keep all steps for a time series on one processor
  chunk[0] = -1;
  chunk[1] = -1;
  chunk[2] = nsteps;
  // Create GA
  p_GA = GA_Create_handle();
  int grp = comm.getGroup();
  GA_Set_data(p_GA,three,dims,C_DBL);
  GA_Set_pgroup(p_GA,grp);
  GA_Set_chunk(p_GA,chunk);
  GA_Allocate(p_GA);
  // Initialize all values to zero
  GA_Zero(p_GA);

  p_comm = comm;
}

/**
 * Basic destructor
 */
gridpack::contingency_analysis::QuantileAnalysis::~QuantileAnalysis()
{
  GA_Destroy(p_GA);
}

/**
 * Save data for a single time step for a single generator
 * @param cfg_idx scenario index for time series
 * @param gen_idx generator index for time series
 * @param vals vector of time series values for a generator
 */
void gridpack::contingency_analysis::QuantileAnalysis::saveData(int cfg_idx,
    int gen_idx, std::vector<double> &vals)
{
  int lo[3], hi[3], ld[2];
  lo[0] = cfg_idx;
  lo[1] = gen_idx;
  lo[2] = 0;
  hi[0] = cfg_idx;
  hi[1] = gen_idx;
  hi[2] = p_nsteps-1;
  ld[0] = 1;
  ld[1] = p_nsteps;
  if (vals.size() > 0) {
    NGA_Put(p_GA,lo,hi,&(vals[0]),ld);
  }

}

/**
 * Save variable names
 * @param name vector of variable names
 */
void gridpack::contingency_analysis::QuantileAnalysis::saveVarNames(
    std::vector<std::string> &names)
{
  p_var_names = names;
  if (names.size() != p_nwatch) {
    printf("Number of variable names does not match number of variables\n");
  }
}

/**
 * Stream data in storage array
 */
void gridpack::contingency_analysis::QuantileAnalysis::writeData()
{
  GA_Sync();
}

/**
 * Calculate quantiles and write them to a file
 * @param quantiles values describing quantiles to be calculated.
 *                  These values should be between 0 and 1.
 * @param dt magnitude time step (in seconds)
 */
void gridpack::contingency_analysis::QuantileAnalysis::exportQuantiles(
    std::vector<double> quantiles, double dt)
{
  int i, j, k, nvals;
  nvals = quantiles.size();
  // Check quantile values to make sure they lie between 0.0 and 1.0
  for (i=0; i<nvals; i++) {
    if (quantiles[i]<0.0 || quantiles[i] > 1.0) {
      //TODO: out of range error
    }
  }
  // Create global array of size n_vals*p_nwatch*p_nsteps to hold nvals
  // quantile values for all watched generators over all timesteps
  int g_quant = GA_Create_handle();
  int dims[3];
  dims[0] = nvals;
  dims[1] = p_nwatch;
  dims[2] = p_nsteps;
  NGA_Set_data(g_quant,3,dims,C_DBL);
  GA_Allocate(g_quant);
  // Sort quantile values so that they run from lowest to highest
  std::sort(quantiles.begin(),quantiles.end());
  std::vector<int> partition(nvals);
  std::vector<double> weight(nvals);

  // find indices that bracket quantiles
  for (i=0; i<nvals; i++) {
    partition[i] = static_cast<int>(quantiles[i]*static_cast<double>(p_nconf));
  }
  if (quantiles[0] == 0.0) partition[0] = 0;
  if (quantiles[nvals-1] == 1.0) partition[nvals-1] = p_nconf-1;
  // calculate weights for each partition
  for (i=0; i<nvals; i++) {
    if (i == 0 && quantiles[0] == 0.0) {
      weight[0] = 1.0;
    } else if (i == nvals-1 && quantiles[i] == 1.0) {
      weight[nvals-1] = 1.0;
    } else {
      double ratio = static_cast<double>(partition[i])/static_cast<double>(p_nconf);
      weight[i] = 1.0 - (ratio - quantiles[i]);
    }
  }
  int lo[3], hi[3], ld[2];
  // nconf, nwatch, nsteps
  std::vector<double> time_slice(p_nconf*p_nwatch);
  lo[0] = 0;
  hi[0] = p_nconf-1;
  lo[1] = 0;
  hi[1] = p_nwatch-1;
  ld[0] = p_nwatch;
  ld[1] = 1;
  // Set up task manager to iterate over steps;
  gridpack::parallel::TaskManager counter(p_comm);
  counter.set(p_nsteps);
  int step;
  int nquar = quantiles.size();
  int rlo[3], rhi[3];
  rlo[0] = 0;
  rlo[1] = 0;
  rhi[0] = nvals-1;
  rhi[1] = p_nwatch-1;
  // loop over steps
  std::vector<double> results(p_nwatch*nvals);
  std::vector<double> qvalues(nvals);
  while (counter.nextTask(&step)){
    lo[2] = step;
    hi[2] = step;
    // Get all data for this step
    NGA_Get(p_GA,lo,hi,&time_slice[0],ld);
    // loop over watched variables
    for (i=0; i<p_nwatch; i++) {
      // get all configuration values for each variable
      std::vector<double> conf(p_nconf);
      for (j=0; j<p_nconf; j++) {
        conf[j] = time_slice[j*p_nwatch+i];
      }
      // Sort values from lowest to highest
      std::sort(conf.begin(),conf.end());
      // Find quantile values 
      for (k=0; k<nvals; k++) {
        int kdx = partition[k];
        if (k == 0 && quantiles[k] == 0.0) {
          qvalues[0] = conf[kdx];
        } else if (k == nvals-1 && quantiles[k] == 1.0) {
          qvalues[k] = conf[kdx];
        } else {
          if (partition[k] < p_nconf-1) {
            qvalues[k] = weight[k]*conf[kdx]+(1.0-weight[k])*conf[kdx+1];
          } else {
            //TODO: report error
          }
        }
      }
      // copy quantile values into results buffer
      for (k=0; k<nvals; k++) {
        results[k*p_nwatch+i] = qvalues[k];
      }
    }
    rlo[2] = step;
    rhi[2] = step;
    // copy results values into global array containing results for all steps
    NGA_Put(g_quant,rlo,rhi,&results[0],ld);
  }
  GA_Sync();
  // Write out results to files. Currently writing out each variable to a
  // separate file
  if (GA_Pgroup_nodeid(p_comm.getGroup()) == 0) {
    // loop over watched variables
    std::vector<double> wbuf(nvals*p_nsteps);
    lo[0] = 0;
    lo[2] = 0;
    hi[0] = nvals-1;
    hi[2] = p_nsteps-1;
    ld[0] = 1;
    ld[1] = p_nsteps;
    for (i=0; i<p_nwatch; i++) {
      lo[1] = i;
      hi[1] = i;
      NGA_Get(g_quant,lo,hi,&wbuf[0],ld);
      FILE *fd = fopen(p_var_names[i].c_str(),"w");
      for (j=0; j<p_nsteps; j++) {
        fprintf(fd,"%16.4f",j*dt);
        for (k=0; k<nvals; k++) {
          fprintf(fd," %16.8f",wbuf[k*p_nsteps+j]);
        }
        fprintf(fd,"\n");
      }
      fclose(fd);
    }
  }
}

/**
 * Basic constructor
 */
gridpack::contingency_analysis::WindDriver::WindDriver(void)
{
}

/**
 * Basic destructor
 */
gridpack::contingency_analysis::WindDriver::~WindDriver(void)
{
}

/**
 * Transfer data from power flow to dynamic simulation
 * @param pf_network power flow network
 * @param ds_network dynamic simulation network
 */
void gridpack::contingency_analysis::WindDriver::transferPFtoDS(
    boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network,
    boost::shared_ptr<DSFullNetwork>
    ds_network)
{
  int numBus = pf_network->numBuses();
  int i;
  gridpack::component::DataCollection *pfData;
  gridpack::component::DataCollection *dsData;
  double rval;
  int ival;
  double temp;
  double pi = 4.0*atan(1.0);
  int itmp = 0;
  for (i=0; i<numBus; i++) {
    pfData = pf_network->getBusData(i).get();
    dsData = ds_network->getBusData(i).get();
    pfData->getValue("BUS_PF_VMAG",&rval);
    dsData->setValue("BUS_PF_VMAG",rval);
    dsData->setValue(BUS_VOLTAGE_MAG,rval);
    temp = rval;
    pfData->getValue("BUS_PF_VANG",&rval);
    dsData->setValue("BUS_PF_VANG",rval);
    dsData->setValue(BUS_VOLTAGE_ANG,rval);
    rval = rval * pi/180.0;
    pfData->getValue("BUS_TYPE", &ival);
    if (ival != 4) {
      itmp ++;
      if ( temp*sin(rval) < 0) {
        printf("Powerflow bus%d mag = %f %fi\n", itmp, temp*cos(rval), temp*sin(rval));
      } else {
        printf("Powerflow bus%d mag = %f +%fi\n", itmp, temp*cos(rval), temp*sin(rval));
      }
    }
    pfData->getValue("LOAD_PL",&rval,0);
    dsData->setValue(LOAD_PL,rval);
    int ngen = 0;
    pfData->getValue(GENERATOR_NUMBER, &ngen);
    if (pfData->getValue(GENERATOR_NUMBER, &ngen)) {
      int j;
      for (j=0; j<ngen; j++) {
        pfData->getValue("GENERATOR_PF_PGEN",&rval,j);
        printf ("Bus:%d, ---PF_PG = %e \n", i, rval);
        dsData->setValue(GENERATOR_PG,rval,j);
        if (ngen >1) printf("save PGEN: %f\n", rval);
        pfData->getValue("GENERATOR_PF_QGEN",&rval,j);
        dsData->setValue(GENERATOR_QG,rval,j);
        if (ngen > 1) printf("save QGEN: %f\n", rval);
      }
    }
  }
}

/**
 * Modify real power parameters in DS network based on values in
 * wind and load files
 * @param genIDs list of bus IDs containing generators to be modified
 * @param genTags list of tags for generators to be modified
 * @param windVals list of new real power values for generators
 * @param loadIDs list of bus IDs containing loads to be modified
 * @param loadTags list of tags for loads to be modified
 * @param loadVals list of new real power values for loads
 * @param pf_network pointer to PF network
 * @param ds_network pointer to DS network
 */
void gridpack::contingency_analysis::WindDriver::resetData(
    std::vector<int> &genIDs, std::vector<std::string> &genTags,
    std::vector<double> &windVals, std::vector<int> &loadIDs,
    std::vector<std::string> &loadTags, std::vector<double> &loadVals,
    boost::shared_ptr<gridpack::powerflow::PFNetwork> &pf_network,
    boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork> &ds_network)
{
  int nsize = genIDs.size();
  int i;
  gridpack::component::DataCollection *data;
  gridpack::dynamic_simulation::DSFullBus *ds_bus;
  gridpack::powerflow::PFBus *pf_bus;
  for (i=0; i<nsize; i++) {
    if (genIDs[i] < 0 || genIDs[i] >= pf_network->numBuses()) {
      printf("Invalid local bus Index: ",genIDs[i]);
    }
    pf_bus = dynamic_cast<gridpack::powerflow::PFBus*>
      (pf_network->getBus(genIDs[i]).get());
    data = pf_network->getBusData(genIDs[i]).get();
    pf_bus->setGeneratorRealPower(genTags[i],windVals[i],data);
    ds_bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
      (ds_network->getBus(genIDs[i]).get());
    data = ds_network->getBusData(genIDs[i]).get();
    ds_bus->setGeneratorRealPower(genTags[i],windVals[i],data);
  }
  
  nsize = loadVals.size();
  for (i=0; i<nsize; i++) {
    if (loadIDs[i] < 0 || loadIDs[i] >= pf_network->numBuses()) {
      printf("Invalid local load Index: ",loadIDs[i]);
    }
    pf_bus = dynamic_cast<gridpack::powerflow::PFBus*>
      (pf_network->getBus(loadIDs[i]).get());
    data = pf_network->getBusData(loadIDs[i]).get();
    pf_bus->setLoadRealPower(loadTags[i],loadVals[i],data);
    ds_bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
      (ds_network->getBus(loadIDs[i]).get());
    data = ds_network->getBusData(loadIDs[i]).get();
    ds_bus->setLoadRealPower(loadTags[i],loadVals[i],data);
  }
}

/**
 * Read faults from external file and form a list of faults
 * @param cursor pointer to open file contain fault or faults
 * @return a list of fault events
 */
std::vector<gridpack::dynamic_simulation::Event>
gridpack::contingency_analysis::WindDriver::
getEvents(gridpack::utility::Configuration::CursorPtr cursor)
{
  gridpack::utility::Configuration::CursorPtr list;
  list = cursor->getCursor("Events");
  gridpack::utility::Configuration::ChildCursors events;
  std::vector<gridpack::dynamic_simulation::Event> ret;
  if (list) {
    list->children(events);
    int size = events.size();
    int idx;
    // Parse fault events
    for (idx=0; idx<size; idx++) {
      gridpack::dynamic_simulation::Event event;
      event.start = events[idx]->get("beginFault",0.0);
      event.end = events[idx]->get("endFault",0.0);
      std::string indices = events[idx]->get("faultBranch","0 0");
      //Parse indices to get from and to indices of branch
      int ntok1 = indices.find_first_not_of(' ',0);
      int ntok2 = indices.find(' ',ntok1);
      if (ntok2 - ntok1 > 0 && ntok1 != std::string::npos && ntok2 !=
          std::string::npos) {
        event.from_idx = atoi(indices.substr(ntok1,ntok2-ntok1).c_str());
        ntok1 = indices.find_first_not_of(' ',ntok2);
        ntok2 = indices.find(' ',ntok1);
        if (ntok1 != std::string::npos && ntok1 < indices.length()) {
          if (ntok2 == std::string::npos) {
            ntok2 = indices.length();
          }
          event.to_idx = atoi(indices.substr(ntok1,ntok2-ntok1).c_str());
        } else {
          event.from_idx = 0;
          event.to_idx = 0;
        }
	event.isBus = false;
	event.isLine = true;
      } else {
        event.from_idx = 0;
        event.to_idx = 0;
      }
      event.step = events[idx]->get("timeStep",0.0);
      if (event.step != 0.0 && event.end != 0.0 &&
          event.from_idx != event.to_idx) {
        ret.push_back(event);
      }
    }
  }
  return ret;
}

/**
 * Read list of branches to monitor in power flow calculation and store
 * contents in internal variables
 * @param cursor pointer to open file containing list of branchs
 */
void gridpack::contingency_analysis::WindDriver::
getWatchLines(gridpack::utility::Configuration::CursorPtr cursor)
{
  gridpack::utility::Configuration::CursorPtr list;
  list = cursor->getCursor("watchBranches");
  gridpack::utility::Configuration::ChildCursors lines;
  p_from.clear();
  p_to.clear();
  p_tags.clear();
  gridpack::utility::StringUtils util;
  if (list) {
    list->children(lines);
    int size = lines.size();
    int idx;
    int from, to;
    std::string tag;
    // Parse watch lines
    for (idx=0; idx<size; idx++) { from = lines[idx]->get("fromBus",-1);
      to = lines[idx]->get("toBus",-1);
      tag = lines[idx]->get("lineTag","-1");
      std::string cleantag = util.clean2Char(tag);
      p_from.push_back(from);
      p_to.push_back(to);
      p_tags.push_back(cleantag);
    }
  }
}

/**
 * Get power flow values from watched branches
 * @param p real power values
 * @param q reactive power values
 * @param network pointer to PF network
 */
void gridpack::contingency_analysis::WindDriver::
getWatchedBranches(std::vector<double> &p, std::vector<double> &q,
    boost::shared_ptr<gridpack::powerflow::PFNetwork> &network)
{
  int nvals = p_from.size();
  int idx, jdx, jvals;
  double *lp, *lq;
  lp = new double[nvals];
  lq = new double[nvals];
  std::vector<int> local;
  gridpack::powerflow::PFBranch *branch;

  // find complex power for transmission elements that are local to this
  // processor
  for (idx=0; idx<nvals; idx++) {
    local = network->getLocalBranchIndices(p_from[idx],p_to[idx]);
    bool found = false;
    jvals = local.size();
    for (jdx=0; jdx<jvals; jdx++) {
      if (network->getActiveBranch(local[jdx])) {
        found = true;
        branch = dynamic_cast<gridpack::powerflow::PFBranch*>(
            network->getBranch(local[jdx]).get());
        gridpack::ComplexType s = branch->getComplexPower(p_tags[idx]);
        lp[idx] = real(s);
        lq[idx] = imag(s);
      }
    }
    if (!found) {
      lp[idx] = 0.0;
      lq[idx] = 0.0;
    }
  }

  // sum up active power over all processors so that every processor has a copy
  // of the data
  network->communicator().sum(lp,nvals);
  network->communicator().sum(lq,nvals);

  // copy data to p and q vectors
  p.clear();
  q.clear();
  for (idx=0; idx<nvals; idx++) {
    p.push_back(lp[idx]);
    q.push_back(lq[idx]);
  }

  delete [] lp;
  delete [] lq;
}

/**
 * Execute application
 */
void gridpack::contingency_analysis::WindDriver::execute(int argc, char** argv)
{
  gridpack::parallel::Communicator world;

  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Total Application");
  timer->start(t_total);

  // read configuration file
  gridpack::utility::Configuration *config
    = gridpack::utility::Configuration::configuration();
  if (argc >= 2 && argv[1] != NULL) {
    char inputfile[256];
    sprintf(inputfile,"%s",argv[1]);
    config->open(inputfile,world);
  } else {
    config->open("input.xml",world);
  }

  // get size of group (communicator) that individual contingency calculations
  // will run on and create task communicator
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Dynamic_simulation");
  int grp_size;
  double Vmin, Vmax;
  if (!cursor->get("groupSize",&grp_size)) {
    grp_size = 1;
  }
  if (grp_size > world.size()) grp_size = world.size();
  if (world.rank() == 0) {
    printf("\nSize of individual task groups: %d\n",grp_size);
  }
  gridpack::parallel::Communicator task_comm = world.divide(grp_size);

  // Read in scenario data for wind generation and loads
  std::string windFile;
  std::string loadFile;
  bool use_loads = true;
  cursor->get("windFile",&windFile);
  use_loads = cursor->get("loadFile",&loadFile);

  cursor = config->getCursor("Configuration.Powerflow");
  bool useNonLinear = false;
  useNonLinear = cursor->get("UseNonLinear", useNonLinear);
  // Find branches to monitor in power flow calculation
  getWatchLines(cursor);
  if (world.rank() == 0) {
    printf("List size %d %d %d\n",p_from.size(),p_to.size(),p_tags.size());
  }
  if (p_from.size()==0 || p_to.size()==0 || p_tags.size()==0) {
    p_watch_lines = false;
  } else {
    p_watch_lines = true;
  }

  // Create global store object to store power flow results
  gridpack::parallel::GlobalStore<double> pf_results(world);
  // local store for powerflow results
  std::vector<double> pf_p;
  std::vector<double> pf_q;

  // set up and run powerflow calculation on sub-communicator. This task is
  // repeated on each sub-communicator
  int t_init_pf = timer->createCategory("Initial Power Flow Simulation");
  timer->start(t_init_pf);
  boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network(new gridpack::powerflow::PFNetwork(task_comm));

  gridpack::powerflow::PFAppModule pf_app;
  pf_app.readNetwork(pf_network, config);
  pf_app.initialize();
  /* TODO: Is this needed? */
#if 1
  if (useNonLinear) {
    pf_app.nl_solve();
  } else {
    pf_app.solve();
  }
  pf_app.write();
  pf_app.saveData();

  getWatchedBranches(pf_p, pf_q, pf_network);
  // this calculation is replicated on all task communicators so only pick
  // results from one replica
  if (world.rank()==0 && p_watch_lines) {
    pf_results.addVector(0,pf_p);
    pf_results.addVector(1,pf_q);
  }
#endif

 
  timer->stop(t_init_pf);
  // Create dynamic simulation applications on each task communicator
  int t_init = timer->createCategory("Initialize Dynamic Simulation");
  timer->start(t_init);
  p_network.reset(new DSFullNetwork(task_comm));
  gridpack::dynamic_simulation::DSFullApp ds_app(task_comm);
  pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
    gridpack::dynamic_simulation::DSFullBranch>(p_network);

  // transfer results from PF calculation to DS calculation
  transferPFtoDS(pf_network, p_network);

  ds_app.setNetwork(p_network,config);
  ds_app.readGenerators();
  ds_app.initialize();
  timer->stop(t_init);
  ds_app.saveTimeSeries(true);
  /* turn and generator watch and set file name */
  ds_app.setGeneratorWatch("watch.txt");
  std::vector<int> bus_ids;
  std::vector<std::string> gen_ids;
  ds_app.getListWatchedGenerators(bus_ids, gen_ids);

  // Read in faults
  cursor = config->getCursor("Configuration.Dynamic_simulation");
  std::string faultfile;
  double total_time, time_step;

  // Find out number of time steps in each simulation
  total_time = cursor->get("simulationTime",0.0);
  time_step = cursor->get("timeStep",0.0);
  if (total_time == 0.0 || time_step == 0.0) {
    // Some kind of error
  } 
  int nsteps = static_cast<int>(total_time/time_step);
  if (world.rank() == 0) {
    printf(" Number of time steps: %d\n",nsteps);
  }

  if (!cursor->get("faultList",&faultfile)) {
    faultfile = "faults.xml";
  }
  // Read in quantile values
  gridpack::utility::StringUtils util;
  std::string quantiles_str;
  if (!cursor->get("quantiles",&quantiles_str)) {
    quantiles_str = "0.0 0.25 0.5 0.75 1.0";
  }
  std::vector<std::string> tokens;
  tokens = util.blankTokenizer(quantiles_str);
  int i, j;
  std::vector<double> quantiles;
  for (i=0; i<tokens.size(); i++) {
    quantiles.push_back(atof(tokens[i].c_str()));
  }
  if (!config->open(faultfile,world) && world.rank() == 0) {
    printf("\nUnable to open fault file: %s\n",faultfile.c_str());
  } else if (world.rank() == 0) {
    printf("\nFaults located in file: %s\n",faultfile.c_str());
  }


  // Find number of generators being watched
  gridpack::utility::Configuration::CursorPtr list;
  list = cursor->getCursor("generatorWatch");
  gridpack::utility::Configuration::ChildCursors watch;
  int num_watch_gen;
  if (list) {
    list->children(watch);
    num_watch_gen = watch.size();
  }

  // get a list of faults
  int t_flts = timer->createCategory("Read Faults");
  timer->start(t_flts);
  cursor = config->getCursor("FaultList.Dynamic_simulation");
  std::vector<gridpack::dynamic_simulation::Event>
    faults = getEvents(cursor);

  if (world.rank() == 0) {
    int idx;
    printf("Number of faults: %d\n",faults.size());
    for (idx = 0; idx < faults.size(); idx++) {
      printf("Fault %d\n",idx);
      printf(" Begin fault: %12.6f End fault: %12.6f\n",
          faults[idx].start, faults[idx].end);
      printf(" From bus: %8d  To bus: %8d\n",
          faults[idx].from_idx,faults[idx].to_idx);
      printf(" Time step: %12.6f\n",faults[idx].step);
    }
  }
  timer->stop(t_flts);
  int me = world.rank();
  /* TODO: Is this needed? */
#if 0
  ds_app.solve(faults[0]);
#endif
 
  // read in wind and load tables
  gridpack::bus_table::BusTable<DSFullNetwork> wind(p_network);
  gridpack::bus_table::BusTable<DSFullNetwork> load(p_network);
  if (!wind.readTable(windFile)) {
    if (world.rank()==0)
      printf("Unable to open wind file: %s\n",windFile.c_str());
  }
  if (use_loads) {
    if (!load.readTable(loadFile)) {
      use_loads = false;
      if (world.rank()==0)
        printf("Unable to open load file: %s\n",loadFile.c_str());
    }
  }
  std::vector<int> busIDs;
  std::vector<int> loadIDs;
  std::vector<std::string> busTags;
  std::vector<std::string> loadTags;
  std::vector<double> windVals;
  std::vector<double> loadVals;
  wind.getLocalIndices(busIDs);
  wind.getTags(busTags);
  if (use_loads) {
    load.getLocalIndices(loadIDs);
    load.getTags(loadTags);
  }

  int numConfigs = wind.getNumColumns();
  if (load.getNumColumns() != wind.getNumColumns()) {
    if (world.rank() == 0) {
      printf("Number of columns in wind file (%d) do not match load file (%d)\n",
          numConfigs, load.getNumColumns());
      printf("Loads are ignored\n");
    }
    use_loads = false;
  }

  
  // set up task manager
  gridpack::parallel::TaskManager taskmgr(world);
  int ntasks = faults.size();
  taskmgr.set(ntasks*numConfigs);

  // Create distributed storage object
  gridpack::contingency_analysis::QuantileAnalysis analysis(world,
      4*bus_ids.size(),ntasks*numConfigs,nsteps);
  // Construct variable names
  std::vector<std::string> var_names;
  char sbuf[128];
  for (i=0; i<gen_ids.size(); i++) {
    std::string tag = gen_ids[i];
    util.trim(tag);
    sprintf(sbuf, "%d_%s_ANG", bus_ids[i], tag.c_str());
    var_names.push_back(sbuf);
    sprintf(sbuf, "%d_%s_SPD", bus_ids[i], tag.c_str());
    var_names.push_back(sbuf);
    sprintf(sbuf, "%d_%s_GENP", bus_ids[i], tag.c_str());
    var_names.push_back(sbuf);
    sprintf(sbuf, "%d_%s_GENQ", bus_ids[i], tag.c_str());
    var_names.push_back(sbuf);
  }
  analysis.saveVarNames(var_names);
  world.barrier();

  // evaluate faults
  int task_id;
  int t_solve = timer->createCategory("Solve Dynamic Simulation");
  int t_file = timer->createCategory("Dynamic Simulation IO");
  while (taskmgr.nextTask(task_comm, &task_id)) {
    int ncnfg = task_id%numConfigs;
    int nfault = (task_id-ncnfg)/numConfigs;
    printf("p[%d] Executing task %d scenario: %d contingency: %d\n",
        world.rank(),task_id,ncnfg,nfault);
    wind.getValues(ncnfg,windVals);
    if (use_loads) {
      load.getValues(ncnfg,loadVals);
    } else {
      loadVals.clear();
    }
    pf_app.reload();
    resetData(busIDs,busTags,windVals,loadIDs,loadTags,loadVals,
	      pf_network,p_network);
    // Recalculate powerflow for new values of generators and loads
 
    if (useNonLinear) {
      pf_app.nl_solve();
    } else {
      pf_app.solve();
    }
    
    pf_app.write();
    pf_app.saveData();

    getWatchedBranches(pf_p, pf_q, pf_network);
    if (task_comm.rank() == 0)
      printf("PF_P size: %d PF_Q size: %d\n",pf_p.size(),pf_q.size());
    // only store results from process 0 on task communicator
    if (task_comm.rank()==0 && p_watch_lines) {
      pf_results.addVector(2*task_id+2,pf_p);
      pf_results.addVector(2*task_id+3,pf_q);
    }

    // transfer results from PF calculation to DS calculation
//    setDSConfig(busIDs,busTags,windVals,loadIDs,loadTags,loadVals,p_network);
    transferPFtoDS(pf_network, p_network);
    ds_app.reload();


    timer->start(t_file);
    sprintf(sbuf,"Event_scn_%d_flt_%d_%d_%d.out",ncnfg,nfault,faults[nfault].from_idx,
        faults[nfault].to_idx);
    ds_app.open(sbuf);
    timer->stop(t_file);
    // Save off time series to watch file
    sprintf(sbuf,"Gen_watch_scn_%d_flt_%d.csv",ncnfg,nfault);
    ds_app.setGeneratorWatch(sbuf);
//    ds_app.solvePreInitialize(faults[nfault]);
    timer->start(t_solve);

    ds_app.solve(faults[nfault]);

    /*
    while(!ds_app.isDynSimuDone()) {
      ds_app.executeOneSimuStep();
    }
    */
    timer->stop(t_solve);
    timer->start(t_file);
    ds_app.write("default");
    std::vector<std::vector<double> > all_series;
    all_series = ds_app.getGeneratorTimeSeries();
    std::vector<int> gen_idx = ds_app.getTimeSeriesMap();
    int iseries;
    for (iseries=0; iseries<gen_idx.size(); iseries++) {
      analysis.saveData(task_id, gen_idx[iseries],all_series[iseries]);
    }
    ds_app.close();
    timer->stop(t_file);
  }
  taskmgr.printStats();

  // Write out powerflow results from head node
  if (p_watch_lines) pf_results.upload();
  if (world.rank()==0 && p_tags.size() > 0) {
    int nvals = p_tags.size()*(ntasks*numConfigs+1);
    double *lp = new double[nvals];
    double *lq = new double[nvals];
    int ncols = ntasks*numConfigs+1;
    int i, j, k;
    int nsize = p_tags.size();
    double *p_ptr = lp;
    double *q_ptr = lq;
    // copy all values from global store to local buffer
    for (i=0; i<ncols; i++) {
      if (p_watch_lines) pf_results.getVector(i,pf_p);
      for (j=0; j<nsize; j++) {
        p_ptr[j] = pf_p[j];
      }
      if (p_watch_lines) pf_results.getVector(i,pf_q);
      for (j=0; j<nsize; j++) {
        q_ptr[j] = pf_q[j];
      }
      p_ptr += nsize;
      q_ptr += nsize;
    }
    // write out values to file
    cursor = config->getCursor("Configuration.Powerflow");
    std::string watchfile;
    cursor->get("watchLineFile",&watchfile);
    printf("Watchfile: (%s)\n",watchfile.c_str());
    if (watchfile.length()==0) watchfile="watchline.csv";
    std::ofstream fout;
    char buf[128];
    fout.open(watchfile.c_str());
    sprintf(buf,"    from       to ckt");
    fout << buf;
    for (j=0; j<ncols; j++) {
      char pbuf[128], qbuf[128];
      if (j == 0) {
        sprintf(pbuf,"  P_base_base");
        sprintf(qbuf,"  Q_base_base");
      } else {
        int ncnfg = (j-1)%numConfigs;
        int nfault = (j-1-ncnfg)/numConfigs;
        sprintf(pbuf,"P_%d_%d",ncnfg,nfault);
        sprintf(qbuf,"Q_%d_%d",ncnfg,nfault);
        sprintf(buf,"  P %4d %4d  Q %4d %4d",ncnfg,nfault,ncnfg,nfault);
      }
      int len = strlen(pbuf);
      for (k=0; k<13-len; k++) {
        fout << " ";
      }
      fout << pbuf;
      len = strlen(qbuf);
      for (k=0; k<13-len; k++) {
        fout << " ";
      }
      fout << qbuf;
    }
    fout << std::endl;
    for (i=0; i<nsize; i++) {
      sprintf(buf,"%8d %8d  %2s",p_from[i],p_to[i],p_tags[i].c_str());
      fout << buf;
      for (j=0; j<ncols; j++) {
        sprintf(buf," %12.6f %12.6f",lp[i+j*nsize],lq[i+j*nsize]);
      fout << buf;
      }
      fout << std::endl;
    }
    fout.close();
    delete [] lp;
    delete [] lq;
  }

  analysis.exportQuantiles(quantiles, time_step);
  timer->stop(t_total);
  if (ntasks*numConfigs >= world.size()/task_comm.size()) {
    timer->dump();
  }
  p_network.reset();
  pf_network.reset();
}
