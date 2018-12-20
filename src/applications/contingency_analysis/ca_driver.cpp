/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ca_driver.cpp
 * @author Bruce Palmer
 * @date   2017-12-08 13:12:46 d3g096
 *
 * @brief Driver for contingency analysis calculation that make use of the
 *        powerflow module to implement individual power flow simulations for
 *        each contingency. The different contingencies are distributed across
 *        separate communicators using the task manager.
 *
 *
 */
// -------------------------------------------------------------

#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "ca_driver.hpp"

#define USE_SUCCESS
#define USE_STATBLOCK
// Sets up multiple communicators so that individual contingency calculations
// can be run concurrently

/**
 * Basic constructor
 */
gridpack::contingency_analysis::CADriver::CADriver(void)
{
}

/**
 * Basic destructor
 */
gridpack::contingency_analysis::CADriver::~CADriver(void)
{
}

/**
 * Get list of contingencies from external file
 * @param cursor pointer to contingencies in input deck
 * @return vector of contingencies
 */
std::vector<gridpack::powerflow::Contingency>
  gridpack::contingency_analysis::CADriver::getContingencies(
      gridpack::utility::Configuration::ChildCursors contingencies)
{
  // The contingencies ChildCursors argument is a vector of configuration
  // pointers. Each element in the vector is pointing at a seperate Contingency
  // block within the Contingencies block in the input file.
  std::vector<gridpack::powerflow::Contingency> ret;
  int size = contingencies.size();
  int i, idx;
  // Create string utilities object to help parse file
  gridpack::utility::StringUtils utils;
  // Loop over all child cursors
  for (idx = 0; idx < size; idx++) {
    std::string ca_type;
    contingencies[idx]->get("contingencyType",&ca_type);
    // Contingency name is used to direct output to different files for each
    // contingency
    std::string ca_name;
    contingencies[idx]->get("contingencyName",&ca_name);
    if (ca_type == "Line") {
      std::string buses;
      contingencies[idx]->get("contingencyLineBuses",&buses);
      std::string names;
      contingencies[idx]->get("contingencyLineNames",&names);
      // Tokenize bus string to get a list of individual buses
      std::vector<std::string> string_vec = utils.blankTokenizer(buses);
      // Convert buses from character strings to ints
      std::vector<int> bus_ids;
      for (i=0; i<string_vec.size(); i++) {
        bus_ids.push_back(atoi(string_vec[i].c_str()));
      }
      string_vec.clear();
      // Tokenize names string to get a list of individual line tags
      string_vec = utils.blankTokenizer(names);
      std::vector<std::string> line_names;
      // clean up line tags so that they are exactly two characters
      for (i=0; i<string_vec.size(); i++) {
        line_names.push_back(utils.clean2Char(string_vec[i]));
      }
      // Check to make sure we found everything
      if (bus_ids.size() == 2*line_names.size()) {
        // Add contingency parameters to contingency struct
        gridpack::powerflow::Contingency contingency;
        contingency.p_name = ca_name;
        contingency.p_type = Branch;
        int i;
        for (i = 0; i < line_names.size(); i++) {
          contingency.p_from.push_back(bus_ids[2*i]);
          contingency.p_to.push_back(bus_ids[2*i+1]);
          contingency.p_ckt.push_back(line_names[i]);
          contingency.p_saveLineStatus.push_back(true);
        }
        // Add branch contingency to contingency list
        ret.push_back(contingency);
      }
    } else if (ca_type == "Generator") {
      std::string buses;
      contingencies[idx]->get("contingencyBuses",&buses);
      std::string gens;
      contingencies[idx]->get("contingencyGenerators",&gens);
      // Tokenize bus string to get a list of individual buses
      std::vector<std::string> string_vec = utils.blankTokenizer(buses);
      std::vector<int> bus_ids;
      // Convert buses from character strings to ints
      for (i=0; i<string_vec.size(); i++) {
        bus_ids.push_back(atoi(string_vec[i].c_str()));
      }
      string_vec.clear();
      // Tokenize gens string to get a list of individual generator tags
      string_vec = utils.blankTokenizer(gens);
      std::vector<std::string> gen_ids;
      // clean up generator tags so that they are exactly two characters
      for (i=0; i<string_vec.size(); i++) {
        gen_ids.push_back(utils.clean2Char(string_vec[i]));
      }
      // Check to make sure we found everything
      if (bus_ids.size() == gen_ids.size()) {
        gridpack::powerflow::Contingency contingency;
        contingency.p_name = ca_name;
        contingency.p_type = Generator;
        int i;
        for (i = 0; i < bus_ids.size(); i++) {
          contingency.p_busid.push_back(bus_ids[i]);
          contingency.p_genid.push_back(gen_ids[i]);
          contingency.p_saveGenStatus.push_back(true);
        }
        // Add generator contingency to contingency list
        ret.push_back(contingency);
      }
    }
  }
  return ret;
}

/**
 * Execute application. argc and argv are standard runtime parameters
 */
void gridpack::contingency_analysis::CADriver::execute(int argc, char** argv)
{
  // Create world communicator for entire simulation
  gridpack::parallel::Communicator world;

  // Get timer instance for timing entire calculation
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Total Application");
  timer->start(t_total);

  // Read configuration file (user specified, otherwise assume that it is
  // call input.xml)
  gridpack::utility::Configuration *config
    = gridpack::utility::Configuration::configuration();
  if (argc >= 2 && argv[1] != NULL) {
    char inputfile[256];
    sprintf(inputfile,"%s",argv[1]);
    config->open(inputfile,world);
  } else {
    config->open("input.xml",world);
  }

  // Get size of group (communicator) that individual contingency calculations
  // will run on and create a task communicator. Each process is part of only
  // one task communicator, even though the world communicator is broken up into
  // many task communicators
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Contingency_analysis");
  int grp_size;
  double Vmin, Vmax;
  // Check to find out if files should be printed for individual power flow
  // calculations
  bool print_calcs;
  std::string tmp_bool;
  gridpack::utility::StringUtils util;
  if (!cursor->get("printCalcFiles",&tmp_bool)) {
    print_calcs = true;
  } else {
    util.toLower(tmp_bool);
    if (tmp_bool == "false") {
      print_calcs = false;
    } else {
      print_calcs = true;
    }
  }
  if (!cursor->get("groupSize",&grp_size)) {
    grp_size = 1;
  }
  if (!cursor->get("minVoltage",&Vmin)) {
    Vmin = 0.9;
  }
  if (!cursor->get("maxVoltage",&Vmax)) {
    Vmax = 1.1;
  }
  // Check for Q limit violations
  bool check_Qlim;
  if (!cursor->get("checkQLimit",&check_Qlim)) {
    check_Qlim = false;
  }
  gridpack::parallel::Communicator task_comm = world.divide(grp_size);

  // Keep track of failed calculations
#ifdef USE_SUCCESS
  std::vector<int> contingency_idx;
  std::vector<bool> contingency_success;
  gridpack::parallel::GlobalVector<bool> ca_success(world);
  std::vector<int> contingency_violation;
  gridpack::parallel::GlobalVector<int> ca_violation(world);
#endif

  // Create powerflow applications on each task communicator
  boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network(new gridpack::powerflow::PFNetwork(task_comm));
  gridpack::powerflow::PFAppModule pf_app;
  // Read in the network from an external file and partition it over the
  // processors in the task communicator. This will read in power flow
  // parameters from the Powerflow block in the input
  pf_app.readNetwork(pf_network,config);
  // Finish initializing the network
  pf_app.initialize();
  //  Set minimum and maximum voltage limits on all buses
  pf_app.setVoltageLimits(Vmin, Vmax);
  // Solve the base power flow calculation. This calculation is replicated on
  // all task communicators
  pf_app.solve();
  // Check for Qlimit violations
  if (check_Qlim && !pf_app.checkQlimViolations()) {
    pf_app.solve();
  }
  // Some buses may violate the voltage limits in the base problem. Flag these
  // buses to ignore voltage violations on them.
  pf_app.ignoreVoltageViolations();

  // Read in contingency file name
  std::string contingencyfile;
  if (!cursor->get("contingencyList",&contingencyfile)) {
    contingencyfile = "contingencies.xml";
  }
  if (world.rank() == 0) printf("Contingency List: %s\n",contingencyfile.c_str());
  // Open contingency file
  bool ok = config->open(contingencyfile,world);

  // Get a list of contingencies. Set cursor so that it points to the
  // Contingencies block in the contingency file
  cursor = config->getCursor(
      "ContingencyList.Contingency_analysis.Contingencies");
  gridpack::utility::Configuration::ChildCursors contingencies;
  if (cursor) cursor->children(contingencies);
  std::vector<gridpack::powerflow::Contingency>
    events = getContingencies(contingencies);
  // Contingencies are now available. Print out a list of contingencies from
  // process 0 (the list is replicated on all processors)
  if (world.rank() == 0) {
    int idx;
    for (idx = 0; idx < events.size(); idx++) {
      printf("Name: %s\n",events[idx].p_name.c_str());
      if (events[idx].p_type == Branch) {
        int nlines = events[idx].p_from.size();
        int j;
        for (j=0; j<nlines; j++) {
          printf(" Line: (from) %d (to) %d (line) \'%s\'\n",
              events[idx].p_from[j],events[idx].p_to[j],
              events[idx].p_ckt[j].c_str());
        }
      } else if (events[idx].p_type == Generator) {
        int nbus = events[idx].p_busid.size();
        int j;
        for (j=0; j<nbus; j++) {
          printf(" Generator: (bus) %d (generator ID) \'%s\'\n",
              events[idx].p_busid[j],events[idx].p_genid[j].c_str());
        }
      }
    }
  }


  // Set up task manager on the world communicator. The number of tasks is
  // equal to the number of contingencies
  gridpack::parallel::TaskManager taskmgr(world);
  int ntasks = events.size();
  taskmgr.set(ntasks);

  int nbus = pf_network->totalBuses();
  // Get bus voltage information for base case
  int i, j;
#ifdef USE_STATBLOCK
  int t_store = timer->createCategory("Store Statistics");
  timer->start(t_store);
  std::vector<std::string> v_vals = pf_app.writeBusString("vr_str");
  int nsize = v_vals.size();
  std::vector<int> mag_ids;
  std::vector<int> ids;
  std::vector<int> branch_ids;
  std::vector<std::string> mag_tags;
  std::vector<std::string> tags;
  std::vector<double> vmag;
  std::vector<double> vang;
  std::vector<int> mag_mask;
  std::vector<int> mask;
  // Find bus IDs and create a dummy tag label and get voltage magnitude
  // and angle for base case
  for (i=0; i<nsize; i++) {
    std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
    int not_isolated = atoi(tokens[3].c_str());
    if (not_isolated == 1) {
      mag_ids.push_back(atoi(tokens[0].c_str()));
      mag_tags.push_back("1 ");
      vmag.push_back(atof(tokens[2].c_str()));
      if (atoi(tokens[4].c_str()) != 0) {
        mag_mask.push_back(2);
      } else {
        mag_mask.push_back(1);
      }
    }
    ids.push_back(atoi(tokens[0].c_str()));
    tags.push_back("1 ");
    vang.push_back(atof(tokens[1].c_str()));
    mask.push_back(1);
  }
  int nmags = vmag.size();
  world.max(&nmags,1);
  world.max(&nbus,1);
#endif
  // Create StatBlock objects for voltage magnitude and angles and add
  // bus IDs to it
#ifdef USE_STATBLOCK
  gridpack::analysis::StatBlock vmag_stats(world,nmags,ntasks+1);
  gridpack::analysis::StatBlock vang_stats(world,nbus,ntasks+1);
#endif
  // Add bus IDs and tags to StatBlock objects as well as base case values of
  // voltage magnitude and angle
#ifdef USE_STATBLOCK
  if (world.rank() == 0) {
    vmag_stats.addRowLabels(mag_ids, mag_tags);
    vang_stats.addRowLabels(ids, tags);
    vmag_stats.addColumnValues(0,vmag,mag_mask);
    vang_stats.addColumnValues(0,vang,mask);
  }
#endif
  // Get generator power information
#ifdef USE_STATBLOCK
  v_vals.clear();
  ids.clear();
  tags.clear();
  mask.clear();
  std::vector<double> pgen;
  std::vector<double> qgen;
  v_vals = pf_app.writeBusString("power");
  nsize = v_vals.size();
  // Find bus IDs and tags for generators and eveluate Pg and Qg for base case
  for (i=0; i<nsize; i++) {
    std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
    if (tokens.size()%4 != 0) {
      printf("Incorrect generator listing\n");
      continue;
    }
    int ngen = tokens.size()/4;
    for (j=0; j<ngen; j++) {
      ids.push_back(atoi(tokens[j*4].c_str()));
      tags.push_back(tokens[j*4+1]);
      pgen.push_back(atof(tokens[j*4+2].c_str()));
      qgen.push_back(atof(tokens[j*4+3].c_str()));
      mask.push_back(1);
    }
  }
  nsize = pgen.size();
  world.max(&nsize,1);
#endif
  // Create StatBlock objects for Pg and Qg and add labels as well as values for
  // base case
#ifdef USE_STATBLOCK
  gridpack::analysis::StatBlock pgen_stats(world,nsize,ntasks+1);
  gridpack::analysis::StatBlock qgen_stats(world,nsize,ntasks+1);
  if (world.rank() == 0) {
    pgen_stats.addRowLabels(ids, tags);
    qgen_stats.addRowLabels(ids, tags);
    pgen_stats.addColumnValues(0,pgen,mask);
    qgen_stats.addColumnValues(0,qgen,mask);
  }
#endif

  // Find flow parameters for all branch lines
#ifdef USE_STATBLOCK
  v_vals.clear();
  ids.clear();
  tags.clear();
  mask.clear();
  std::vector<int> id1;
  std::vector<int> id2;
  std::vector<double> pmin, pmax;
  std::vector<double> pflow;
  std::vector<double> qflow;
  std::vector<double> perf;
  v_vals = pf_app.writeBranchString("flow_str");
  nsize = v_vals.size();
  // Parse branch line endpoints as well as line IDs and values of P and Q for
  // base case
  for (i=0; i<nsize; i++) {
    std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
    if (tokens.size()%8 != 0) {
      printf("Incorrect branch power flow listing\n");
      continue;
    }
    int nline = tokens.size()/8;
    for (j=0; j<nline; j++) {
      id1.push_back(atoi(tokens[j*8].c_str()));
      id2.push_back(atoi(tokens[j*8+1].c_str()));
      tags.push_back(tokens[j*8+2]);
      pflow.push_back(atof(tokens[j*8+3].c_str()));
      qflow.push_back(atof(tokens[j*8+4].c_str()));
      perf.push_back(atof(tokens[j*8+5].c_str()));
      pmin.push_back(-atof(tokens[j*8+6].c_str()));
      pmax.push_back(atof(tokens[j*8+6].c_str()));
      if (atoi(tokens[j*8+7].c_str()) == 0) {
        mask.push_back(1);
      } else {
        mask.push_back(2);
      }
    }
  }
  nsize = pflow.size();
  world.max(&nsize,1);
#endif
  // Create StatBlock objects for flow parameters and add labels and base case
  // values
#ifdef USE_STATBLOCK
  gridpack::analysis::StatBlock pflow_stats(world,nsize,ntasks+1);
  gridpack::analysis::StatBlock qflow_stats(world,nsize,ntasks+1);
  gridpack::analysis::StatBlock perf_stats(world,nsize,ntasks+1);
  if (world.rank() == 0) {
    pflow_stats.addRowLabels(id1, id2, tags);
    qflow_stats.addRowLabels(id1, id2, tags);
    perf_stats.addRowLabels(id1, id2, tags);
    pflow_stats.addColumnValues(0,pflow,mask);
    qflow_stats.addColumnValues(0,qflow,mask);
    perf_stats.addColumnValues(0,perf,mask);
    pflow_stats.addRowMinValue(pmin);
    qflow_stats.addRowMinValue(pmin);
    pflow_stats.addRowMaxValue(pmax);
    qflow_stats.addRowMaxValue(pmax);
  }
  timer->stop(t_store);
#endif
  if (check_Qlim) pf_app.clearQlimViolations();


  // Evaluate contingencies using the task manager
  int task_id;
  char sbuf[128];
  // nextTask returns the same task_id on all processors in task_comm. When the
  // calculation runs out of task, nextTask will return false.
  while (taskmgr.nextTask(task_comm, &task_id)) {
    printf("Executing task %d on process %d\n",task_id,world.rank());
    sprintf(sbuf,"%s.out",events[task_id].p_name.c_str());
    // Open a new file, based on the contingency name, to store results from
    // this particular contingency calculation
    if (print_calcs) pf_app.open(sbuf);
    // Write out information to the top of the output file providing some
    // information on the contingency
    sprintf(sbuf,"\nRunning task on %d processes\n",task_comm.size());
    if (print_calcs) pf_app.writeHeader(sbuf);
    if (events[task_id].p_type == Branch) {
      int nlines = events[task_id].p_from.size();
      int j;
      for (j=0; j<nlines; j++) {
        sprintf(sbuf," Line: (from) %d (to) %d (line) \'%s\'\n",
            events[task_id].p_from[j],events[task_id].p_to[j],
            events[task_id].p_ckt[j].c_str());
        printf("p[%d] Line: (from) %d (to) %d (line) \'%s\'\n",
            pf_network->communicator().rank(),
            events[task_id].p_from[j],events[task_id].p_to[j],
            events[task_id].p_ckt[j].c_str());
      }
    } else if (events[task_id].p_type == Generator) {
      int nbus = events[task_id].p_busid.size();
      int j;
      for (j=0; j<nbus; j++) {
        sprintf(sbuf," Generator: (bus) %d (generator ID) \'%s\'\n",
            events[task_id].p_busid[j],events[task_id].p_genid[j].c_str());
        printf("p[%d] Generator: (bus) %d (generator ID) \'%s\'\n",
            pf_network->communicator().rank(),
            events[task_id].p_busid[j],events[task_id].p_genid[j].c_str());
      }
    }
    if (print_calcs) pf_app.writeHeader(sbuf);
    // Reset all voltages back to their original values
    pf_app.resetVoltages();
    // Set contingency
    pf_app.setContingency(events[task_id]);
    // Solve power flow equations for this system
#ifdef USE_SUCCESS
    contingency_idx.push_back(task_id);
#endif
    if (pf_app.solve()) {
#ifdef USE_SUCCESS
      contingency_success.push_back(true);
#endif
      if (check_Qlim && !pf_app.checkQlimViolations()) {
        pf_app.solve();
      }
      // If power flow solution is successful, write out voltages and currents
      if (print_calcs) pf_app.write();
      // Check for violations
      bool ok1 = pf_app.checkVoltageViolations();
      bool ok2 = pf_app.checkLineOverloadViolations();
      bool ok = ok1 && ok2;
      // Include results of violation checks in output
      if (ok) {
        sprintf(sbuf,"\nNo violation for contingency %s\n",
            events[task_id].p_name.c_str());
#ifdef USE_SUCCESS
        contingency_violation.push_back(1);
#endif
      } 
      if (!ok1) {
        sprintf(sbuf,"\nBus Violation for contingency %s\n",
            events[task_id].p_name.c_str());
      }
      if (print_calcs) pf_app.print(sbuf);
      if (print_calcs) pf_app.writeCABus();
      if (!ok2) {
        sprintf(sbuf,"\nBranch Violation for contingency %s\n",
            events[task_id].p_name.c_str());
      }

#ifdef USE_SUCCESS
      if (!ok1 && !ok2) {
        contingency_violation.push_back(4);
      } else if (!ok1) {
        contingency_violation.push_back(2);
      } else if (!ok2) {
        contingency_violation.push_back(3);
      }
#endif
        
      if (print_calcs) pf_app.print(sbuf);
      if (print_calcs) pf_app.writeCABranch();
      // Get strings of data from power flow calculation and parse them to
      // extract numerical values. Store these values in vectors and then
      // add them to StatBlock objects
#ifdef USE_STATBLOCK
      timer->start(t_store);
      vmag.clear();
      vang.clear();
      mask.clear();
      mag_mask.clear();
      v_vals.clear();
      v_vals = pf_app.writeBusString("vr_str");
      nsize = v_vals.size();
      for (i=0; i<nsize; i++) {
        std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
        int not_isolated = atoi(tokens[3].c_str());
        if (not_isolated == 1) {
          vmag.push_back(atof(tokens[2].c_str()));
          if (atoi(tokens[4].c_str()) != 0) {
            mag_mask.push_back(2);
          } else {
            mag_mask.push_back(1);
          }
        }
        vang.push_back(atof(tokens[1].c_str()));
        mask.push_back(1);
      }
#endif
#ifdef USE_STATBLOCK
      if (task_comm.rank() == 0) {
        vmag_stats.addColumnValues(task_id+1,vmag,mag_mask);
        vang_stats.addColumnValues(task_id+1,vang,mask);
      }
#endif
#ifdef USE_STATBLOCK
      pgen.clear();
      qgen.clear();
      mask.clear();
      v_vals.clear();
      v_vals = pf_app.writeBusString("power");
      nsize = v_vals.size();
      for (i=0; i<nsize; i++) {
        std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
        if (tokens.size()%4 != 0) {
          printf("Incorrect generator listing\n");
          continue;
        }
        int ngen = tokens.size()/4;
        for (j=0; j<ngen; j++) {
          pgen.push_back(atof(tokens[j*4+2].c_str()));
          qgen.push_back(atof(tokens[j*4+3].c_str()));
          mask.push_back(1);
        }
      }
#endif
#ifdef USE_STATBLOCK
      if (task_comm.rank() == 0) {
        pgen_stats.addColumnValues(task_id+1,pgen,mask);
        qgen_stats.addColumnValues(task_id+1,qgen,mask);
      }
#endif
#ifdef USE_STATBLOCK
      pflow.clear();
      qflow.clear();
      perf.clear();
      mask.clear();
      v_vals.clear();
      v_vals = pf_app.writeBranchString("flow_str");
      nsize = v_vals.size();
      for (i=0; i<nsize; i++) {
        std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
        if (tokens.size()%8 != 0) {
          printf("Incorrect branch power flow listing\n");
          continue;
        }
        int nline = tokens.size()/8;
        for (j=0; j<nline; j++) {
          pflow.push_back(atof(tokens[j*8+3].c_str()));
          qflow.push_back(atof(tokens[j*8+4].c_str()));
          perf.push_back(atof(tokens[j*8+5].c_str()));
          if (atoi(tokens[j*8+7].c_str()) == 0) {
            mask.push_back(1);
          } else {
            mask.push_back(2);
          }
        }
      }
#endif
#ifdef USE_STATBLOCK
      if (task_comm.rank() == 0) {
        pflow_stats.addColumnValues(task_id+1,pflow,mask);
        qflow_stats.addColumnValues(task_id+1,qflow,mask);
        perf_stats.addColumnValues(task_id+1,perf,mask);
      }
      timer->stop(t_store);
#endif
      if (check_Qlim) pf_app.clearQlimViolations();
    } else {
#ifdef USE_SUCCESS
      contingency_success.push_back(false);
      contingency_violation.push_back(0);
#endif
      sprintf(sbuf,"\nDivergent for contingency %s\n",
          events[task_id].p_name.c_str());
      if (print_calcs) pf_app.print(sbuf);
      // Add dummy values to StatBlock object. Mask value is set to 0 for all
      // network elements to indicate calculation failure
#ifdef USE_STATBLOCK
      timer->start(t_store);
      vmag.clear();
      vang.clear();
      mask.clear();
      mag_mask.clear();
      v_vals.clear();
      v_vals = pf_app.writeBusString("vfail_str");
      nsize = v_vals.size();
      for (i=0; i<nsize; i++) {
        std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
        int not_isolated = atoi(tokens[3].c_str());
        if (not_isolated == 1) {
          vmag.push_back(0.0);
          mag_mask.push_back(0);
        }
        vang.push_back(0.0);
        mask.push_back(0);
      }
#endif
#ifdef USE_STATBLOCK
      if (task_comm.rank() == 0) {
        vmag_stats.addColumnValues(task_id+1,vmag,mag_mask);
        vang_stats.addColumnValues(task_id+1,vang,mask);
      }
#endif
#ifdef USE_STATBLOCK
      pgen.clear();
      qgen.clear();
      mask.clear();
      v_vals.clear();
      v_vals = pf_app.writeBusString("pfail_str");
      nsize = v_vals.size();
      for (i=0; i<nsize; i++) {
        std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
        if (tokens.size()%4 != 0) {
          printf("Incorrect generator listing\n");
          continue;
        }
        int ngen = tokens.size()/4;
        for (j=0; j<ngen; j++) {
          pgen.push_back(0.0);
          qgen.push_back(0.0);
          mask.push_back(0);
        }
      }
#endif
#ifdef USE_STATBLOCK
      if (task_comm.rank() == 0) {
        pgen_stats.addColumnValues(task_id+1,pgen,mask);
        qgen_stats.addColumnValues(task_id+1,qgen,mask);
      }
#endif
#ifdef USE_STATBLOCK
      pflow.clear();
      qflow.clear();
      perf.clear();
      mask.clear();
      v_vals.clear();
      v_vals = pf_app.writeBranchString("fail_str");
      nsize = v_vals.size();
      for (i=0; i<nsize; i++) {
        std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
        if (tokens.size()%8 != 0) {
          printf("Incorrect branch power flow listing\n");
          continue;
        }
        int nline = tokens.size()/8;
        for (j=0; j<nline; j++) {
          pflow.push_back(0.0);
          qflow.push_back(0.0);
          perf.push_back(0.0);
          mask.push_back(0);
        }
      }
#endif
#ifdef USE_STATBLOCK
      if (task_comm.rank() == 0) {
        pflow_stats.addColumnValues(task_id+1,pflow,mask);
        qflow_stats.addColumnValues(task_id+1,qflow,mask);
        perf_stats.addColumnValues(task_id+1,perf,mask);
      }
      timer->stop(t_store);
#endif
    } 
    // Return network to its original base case state
    pf_app.unSetContingency(events[task_id]);
    // Close output file for this contingency
    if (print_calcs) pf_app.close();
  }
  // Print statistics from task manager describing the number of tasks performed
  // per processor
  taskmgr.printStats();

  // Gather stats on successful contingency calculations
#ifdef USE_SUCCESS
  if (task_comm.rank() == 0) {
    ca_success.addElements(contingency_idx, contingency_success);
    ca_violation.addElements(contingency_idx, contingency_violation);
  }
  ca_success.upload();
  ca_violation.upload();
  // Write out stats on successful calculations
  if (world.rank() == 0) {
    contingency_idx.clear();
    contingency_success.clear();
    contingency_violation.clear();
    for (i=0; i<ntasks; i++) contingency_idx.push_back(i);
    ca_success.getData(contingency_idx, contingency_success);
    contingency_success.clear();
    ca_violation.getData(contingency_idx, contingency_violation);
    std::ofstream fout;
    fout.open("success.txt");
    for (i=0; i<ntasks; i++) {
      if (contingency_success[i]) {
        fout << "contingency: " << i+1 << " success: true";
        if (contingency_violation[i] == 1) {
          fout << " violation: none" << std::endl;
        } else if (contingency_violation[i] == 2) {
          fout << " violation: bus" << std::endl;
        } else if (contingency_violation[i] == 3) {
          fout << " violation: branch" << std::endl;
        } else if (contingency_violation[i] == 4) {
          fout << " violation: bus and branch" << std::endl;
        }
      } else {
        fout << "contingency: " << i+1 << " success: false" << std::endl;
      }
    }
    fout.close();
  }
#endif

  // Print out statistics on contingencies
#ifdef USE_STATBLOCK
  int t_stats = timer->createCategory("Write Statistics");
  timer->start(t_stats);
  vmag_stats.writeMeanAndRMS("vmag.txt",1,false);
  vmag_stats.writeMinAndMax("vmag_mm.txt",1,false);
  if (check_Qlim) vmag_stats.writeMaskValueCount("pq_change_cnt.txt",2,false);
  vang_stats.writeMeanAndRMS("vang.txt",1,false);
  vang_stats.writeMinAndMax("vang_mm.txt",1,false);
  pgen_stats.writeMeanAndRMS("pgen.txt",1);
  pgen_stats.writeMinAndMax("pgen_mm.txt",1);
  qgen_stats.writeMeanAndRMS("qgen.txt",1);
  qgen_stats.writeMinAndMax("qgen_mm.txt",1);
  pflow_stats.writeMeanAndRMS("pflow.txt",1);
  pflow_stats.writeMinAndMax("pflow_mm.txt",1);
  pflow_stats.writeMaskValueCount("line_flt_cnt.txt",2);
  qflow_stats.writeMeanAndRMS("qflow.txt",1);
  qflow_stats.writeMinAndMax("qflow_mm.txt",1);
  perf_stats.writeMinAndMax("perf_mm.txt",1);
  perf_stats.sumColumnValues("perf_sum.txt",1);
  timer->stop(t_stats);
#endif
  timer->stop(t_total);
  // If all processors executed at least one task, then print out timing
  // statistics (this printout does not work if some processors do not define
  // all timing variables)
  if (contingencies.size()*grp_size >= world.size()) {
    timer->dump();
  }
}

