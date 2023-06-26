/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wind_driver2.cpp
 *
 * @brief Alternate version of wind driver that can be used to
 * assess different types of events
 *
 *
 */
// -------------------------------------------------------------

#include <algorithm>
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/dynamic_simulation_full_y/dsf_app_module.hpp"
#include "wind_driver.hpp"

/**
 * Execute application
 */
void gridpack::contingency_analysis::WindDriver::execute2(int argc, char** argv)
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

  bool      pf_converged;
  gridpack::powerflow::PFAppModule pf_app;
  pf_app.readNetwork(pf_network, config);
  pf_app.initialize();

  if (useNonLinear) {
    pf_app.nl_solve();
  } else {
    pf_app.solve();
  }
  //  pf_app.write();
  pf_app.saveData();

  getWatchedBranches(pf_p, pf_q, pf_network);
  // this calculation is replicated on all task communicators so only pick
  // results from one replica
  if (world.rank()==0 && p_watch_lines) {
    pf_results.addVector(0,pf_p);
    pf_results.addVector(1,pf_q);
  }

  timer->stop(t_init_pf);
  // Create dynamic simulation applications on each task communicator
  int t_init = timer->createCategory("Initialize Dynamic Simulation");
  timer->start(t_init);
  p_network.reset(new DSFullNetwork(task_comm));
  gridpack::dynamic_simulation::DSFullApp ds_app(task_comm);
  
  pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
    gridpack::dynamic_simulation::DSFullBranch>(p_network);

  // transfer results from PF calculation to DS calculation
  ds_app.transferPFtoDS(pf_network, p_network);

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

  // Find out number of time steps in each simulation
  double total_time, time_step;
  total_time = ds_app.getFinalTime();
  time_step = ds_app.getTimeStep();
  if (total_time == 0.0 || time_step == 0.0) {
    // Some kind of error
  } 
  int nsteps = static_cast<int>(total_time/time_step);
  if (world.rank() == 0) {
    printf(" Number of time steps: %d\n",nsteps);
  }

  cursor = config->getCursor("Configuration.Dynamic_simulation");

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

    /* Get list of generator watch */
  std::string watchlistfile;
  gridpack::utility::Configuration *config_watch
    = gridpack::utility::Configuration::configuration();

  gridpack::utility::Configuration::CursorPtr cursor_watch;
  if (cursor->get("generatorWatchList",&watchlistfile)) {
    if (!config_watch->open(watchlistfile,world) && world.rank() == 0) {
      printf("\nUnable to open watchlist file: %s\n",watchlistfile.c_str());
    } else if (world.rank() == 0) {
      printf("\nWatchlist located in file: %s\n",watchlistfile.c_str());
    }

    cursor_watch = config_watch->getCursor("Dynamic_simulation.generatorWatch");
    ds_app.setGeneratorWatch(cursor_watch);
    ds_app.getListWatchedGenerators(bus_ids, gen_ids);
  } else {
    // Find number of generators being watched
    gridpack::utility::Configuration::CursorPtr list;
    list = cursor->getCursor("generatorWatch");
    gridpack::utility::Configuration::ChildCursors watch;
    int num_watch_gen;
    if (list) {
      list->children(watch);
      num_watch_gen = watch.size();
    }
  }

  // Read in events
  std::string eventfile;

  if (!cursor->get("EventList",&eventfile)) {
    eventfile = "faults.xml";
  }

  if (!config->open(eventfile,world) && world.rank() == 0) {
    printf("\nUnable to open event file: %s\n",eventfile.c_str());
  } else if (world.rank() == 0) {
    printf("\nEvent located in file: %s\n",eventfile.c_str());
  }

  // get a list of faults
  int t_evts = timer->createCategory("Read Events");
  timer->start(t_evts);
  cursor = config->getCursor("EventList.Dynamic_simulation");
  std::vector<gridpack::dynamic_simulation::Event>
    events = ds_app.getEvents(cursor);

  if (world.rank() == 0) {
    int idx;
    printf("Number of events: %d\n",events.size());
    for (idx = 0; idx < events.size(); idx++) {
      printf("Fault %d\n",idx);
      printf(" Begin fault: %12.6f End fault: %12.6f\n",
          events[idx].start, events[idx].end);
      if (events[idx].isBusFault) {
        printf(" Bus ID: %8d\n", events[idx].bus_idx);
      } else if (events[idx].isLineStatus) {
        printf(" From bus: %d To bus: %d\n",
            events[idx].from_idx, events[idx].to_idx);
      } else if (events[idx].isGenStatus) {
      printf(" Gen ID: %s Bus ID: %d\n",
          events[idx].tag.c_str(),events[idx].bus_idx);
      }
    }
  }
  timer->stop(t_evts);
  
  int me = world.rank();
 
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
  if (use_loads) {
    if (load.getNumColumns() != wind.getNumColumns()) {
      if (world.rank() == 0) {
        printf("Number of columns in wind file (%d) do not match load file (%d)\n",
            numConfigs, load.getNumColumns());
        printf("Loads are ignored\n");
      }
      use_loads = false;
    }
  }

  
  // set up task manager
  gridpack::parallel::TaskManager taskmgr(world);
  int ntasks = events.size();
  taskmgr.set(ntasks*numConfigs);

  // Create distributed storage object
  gridpack::contingency_analysis::QuantileAnalysis analysis(world,
      4*bus_ids.size(),ntasks*numConfigs,nsteps-1);
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

  // evaluate events
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
      pf_converged = pf_app.nl_solve();
    } else {
      pf_converged = pf_app.solve();
    }

    if(!pf_converged) {
      printf("p[%d] Power flow diverged %d scenario: %d contingency: %d, skipping this task\n",
	     world.rank(),task_id,ncnfg,nfault);
      continue;
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
    ds_app.transferPFtoDS(pf_network, p_network);
    ds_app.reset();

    timer->start(t_file);
    if (events[nfault].isBusFault) {
      sprintf(sbuf,"Event_scn_%d_flt_%d_bus_%d.out",ncnfg,nfault,
          events[nfault].bus_idx);
    } else if (events[nfault].isLineStatus) {
      sprintf(sbuf,"Event_scn_%d_flt_%d_line_%d_%d.out",ncnfg,nfault,
          events[nfault].from_idx,
          events[nfault].to_idx);
    } else if (events[nfault].isGenStatus) {
      std::string chk = events[nfault].tag;
      util.trim(chk);
      sprintf(sbuf,"Event_scn_%d_flt_%d_gen_%d_%s.out",ncnfg,nfault,
          events[nfault].bus_idx,chk.c_str());
    }
    ds_app.open(sbuf);
    timer->stop(t_file);
    // Save off time series to watch file
    sprintf(sbuf,"Gen_watch_scn_%d_flt_%d.csv",ncnfg,nfault);
    ds_app.setGeneratorWatch(sbuf,cursor_watch);
    timer->start(t_solve);

    ds_app.setup();
    ds_app.setEvent(events[nfault]);
    ds_app.run();

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
