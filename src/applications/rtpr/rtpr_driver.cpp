/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   rtpr_driver.cpp
 * @author Bruce Palmer
 * @date   2019-10-09 13:12:46 d3g293
 *
 * @brief Driver for real-time path rating calculation that make use of the
 *        powerflow module to implement individual power flow simulations for
 *        each contingency. The different contingencies are distributed across
 *        separate communicators using the task manager.
 *
 *
 */
// -------------------------------------------------------------

#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "rtpr_driver.hpp"

#define USE_SUCCESS
// Sets up multiple communicators so that individual contingency calculations
// can be run concurrently

/**
 * Basic constructor
 */
gridpack::rtpr::RTPRDriver::RTPRDriver(void)
{
}

/**
 * Basic destructor
 */
gridpack::rtpr::RTPRDriver::~RTPRDriver(void)
{
}

/**
 * Get list of contingencies from external file
 * @param cursor pointer to contingencies in input deck
 * @return vector of contingencies
 */
std::vector<gridpack::powerflow::Contingency>
  gridpack::rtpr::RTPRDriver::getContingencies(
      gridpack::utility::Configuration::ChildCursors &contingencies)
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
 * Dummy struct to get around some problems with vectors of characters
 */
struct char2 {
  char str[3];
};

/**
 * Create a list of all N-1 generator contingencies for a given area
 * @param network power grid network on which contingencies are defined 
 * @param area index of area that will generate contingencies
 * @return vector of contingencies
 */
std::vector<gridpack::powerflow::Contingency> 
  gridpack::rtpr::RTPRDriver::createGeneratorContingencies(
  boost::shared_ptr<gridpack::powerflow::PFNetwork> network, int area)
{
  std::vector<gridpack::powerflow::Contingency> ret;
  gridpack::utility::StringUtils util;
  int nbus = network->numBuses();
  int i,j,iarea;
  std::vector<int> bus_ids;
  std::vector<char2> tags;
  char2 buf;
  // Loop over all buses
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      // Get data collection object
      gridpack::component::DataCollection *data
        = network->getBusData(i).get();
      int ngen=0; 
      data->getValue(BUS_AREA, &iarea);
      if (iarea == area) {
        if (data->getValue(GENERATOR_NUMBER, &ngen)) {
          for (j=0; j<ngen; j++) {
            // Check generator status
            int status = 0;
            data->getValue(GENERATOR_STAT, &status, j);
            if (status != 0) {
              // Generator is active in base case. Add a contingency
              bus_ids.push_back(network->getOriginalBusIndex(i));
              std::string tag, clean_tag;
              data->getValue(GENERATOR_ID,&tag,j);
              clean_tag = util.clean2Char(tag);
              strncpy(buf.str,clean_tag.c_str(),2);
              buf.str[2] = '\0';
              tags.push_back(buf);
            }
          }
        }
      }
    }
  }
  // Have a list of local faults. Find out total number of faults.
  int nflt = bus_ids.size();
  int nproc = network->communicator().size();
  int me = network->communicator().rank();
  std::vector<int> sizes(nproc);
  for (i=0; i<nproc; i++) sizes[i] = 0;
  sizes[me] = nflt;
  network->communicator().sum(&sizes[0], nproc);
  // create global list containing buses and generator tags
  nflt = 0;
  for (i=0; i<nproc; i++) nflt += sizes[i];
  int g_bus = GA_Create_handle();
  int one = 1;
  NGA_Set_data(g_bus,one,&nflt,C_INT);
  NGA_Allocate(g_bus);
  int g_ids = GA_Create_handle();
  int idType = NGA_Register_type(sizeof(char2));
  NGA_Set_data(g_ids,1,&nflt,idType);
  NGA_Allocate(g_ids);
  int lo, hi;
  lo = 0; 
  for (i=0; i<me; i++) lo += sizes[i];
  hi = lo + sizes[me] - 1;
  // Create a complete list on all processes
  NGA_Put(g_bus,&lo,&hi,&bus_ids[0],&one);
  NGA_Put(g_ids,&lo,&hi,&tags[0],&one);
  network->communicator().sync();
  bus_ids.resize(nflt);
  tags.resize(nflt);
  lo = 0;
  hi = nflt-1;
  NGA_Get(g_bus,&lo,&hi,&bus_ids[0],&one);
  NGA_Get(g_ids,&lo,&hi,&tags[0],&one);
  GA_Destroy(g_bus);
  GA_Destroy(g_ids);
  NGA_Deregister_type(idType);
  // create a list of contingencies
  char sbuf[128];
  for (i=0; i<nflt; i++) {
    gridpack::powerflow::Contingency fault;
    fault.p_type = Generator;
    fault.p_busid.push_back(bus_ids[i]);
    std::string tag = tags[i].str;
    fault.p_genid.push_back(util.clean2Char(tag));
    fault.p_saveGenStatus.push_back(true);
    sprintf(sbuf,"GenCTG%d",i+1);
    fault.p_name = sbuf;
    ret.push_back(fault);
    if (GA_Nodeid() == 0) {
      printf("bus: %d tag: (%s) name: (%s)\n",fault.p_busid[0],
          fault.p_genid[0].c_str(),fault.p_name.c_str());
    }
  }
  return ret;
}

  /**
   * Create a list of all N-1 branch contingencies for a given area
   * @param network power grid network on which contingencies are defined 
   * @param area index of area that will generate contingencies
   * @return vector of contingencies
   */
std::vector<gridpack::powerflow::Contingency>
  gridpack::rtpr::RTPRDriver::createBranchContingencies(
  boost::shared_ptr<gridpack::powerflow::PFNetwork> network, int area)
{
  std::vector<gridpack::powerflow::Contingency> ret;
  gridpack::utility::StringUtils util;
  int nbranch = network->numBranches();
  int i,j,idx1,idx2,i1,i2,area1,area2;
  std::vector<int> bus_from;
  std::vector<int> bus_to;
  std::vector<char2> tags;
  char2 buf;
  // Loop over all branches
  for (i=0; i<nbranch; i++) {
    if (network->getActiveBranch(i)) {
      // Get data collection object
      gridpack::component::DataCollection *data
        = network->getBranchData(i).get();
      int nline=0; 
      network->getBranchEndpoints(i,&i1,&i2);
      network->getBusData(i1)->getValue(BUS_AREA,&area1);
      network->getBusData(i2)->getValue(BUS_AREA,&area2);
      if (area == area1 && area == area2) {
        if (data->getValue(BRANCH_NUM_ELEMENTS, &nline)) {
          for (j=0; j<nline; j++) {
            // Check generator status
            int status = 0;
            data->getValue(BRANCH_STATUS, &status, j);
            if (status != 0) {
              // Generator is active in base case. Add a contingency
              network->getOriginalBranchEndpoints(i,&idx1,&idx2);
              bus_from.push_back(idx1);
              bus_to.push_back(idx2);
              std::string tag, clean_tag;
              data->getValue(BRANCH_CKT,&tag,j);
              clean_tag = util.clean2Char(tag);
              strncpy(buf.str,clean_tag.c_str(),2);
              buf.str[2] = '\0';
              tags.push_back(buf);
            }
          }
        }
      }
    }
  }
  // Have a list of local faults. Find out total number of faults.
  int nflt = bus_from.size();
  int nproc = network->communicator().size();
  int me = network->communicator().rank();
  std::vector<int> sizes(nproc);
  for (i=0; i<nproc; i++) sizes[i] = 0;
  sizes[me] = nflt;
  network->communicator().sum(&sizes[0], nproc);
  // create global list containing buses and generator tags
  nflt = 0;
  for (i=0; i<nproc; i++) nflt += sizes[i];
  int g_from = GA_Create_handle();
  int one = 1;
  NGA_Set_data(g_from,one,&nflt,C_INT);
  NGA_Allocate(g_from);
  int g_to = GA_Create_handle();
  NGA_Set_data(g_to,one,&nflt,C_INT);
  NGA_Allocate(g_to);
  int g_ids = GA_Create_handle();
  int idType = NGA_Register_type(sizeof(char2));
  NGA_Set_data(g_ids,1,&nflt,idType);
  NGA_Allocate(g_ids);
  int lo, hi;
  lo = 0; 
  for (i=0; i<me; i++) lo += sizes[i];
  hi = lo + sizes[me] - 1;
  // Create a complete list on all processes
  NGA_Put(g_from,&lo,&hi,&bus_from[0],&one);
  NGA_Put(g_to,&lo,&hi,&bus_to[0],&one);
  NGA_Put(g_ids,&lo,&hi,&tags[0],&one);
  network->communicator().sync();
  bus_from.resize(nflt);
  bus_to.resize(nflt);
  tags.resize(nflt);
  lo = 0;
  hi = nflt-1;
  NGA_Get(g_from,&lo,&hi,&bus_from[0],&one);
  NGA_Get(g_to,&lo,&hi,&bus_to[0],&one);
  NGA_Get(g_ids,&lo,&hi,&tags[0],&one);
  GA_Destroy(g_from);
  GA_Destroy(g_to);
  GA_Destroy(g_ids);
  NGA_Deregister_type(idType);
  // create a list of contingencies
  char sbuf[128];
  for (i=0; i<nflt; i++) {
    gridpack::powerflow::Contingency fault;
    fault.p_type = Branch;
    fault.p_from.push_back(bus_from[i]);
    fault.p_to.push_back(bus_to[i]);
    std::string tag = tags[i].str;
    fault.p_ckt.push_back(util.clean2Char(tag));
    fault.p_saveLineStatus.push_back(true);
    sprintf(sbuf,"LineCTG%d",i+1);
    fault.p_name = sbuf;
    ret.push_back(fault);
    if (GA_Nodeid() == 0) {
      printf("from: %d to: %d tag: (%s) name: (%s)\n",fault.p_from[0],
          fault.p_to[0],
          fault.p_ckt[0].c_str(),fault.p_name.c_str());
    }
  }
  return ret;
}

/**
 * Get list of tie lines
 * @param cursor pointer to tie lines in input deck
 * @return vector of contingencies
 */
std::vector<gridpack::rtpr::TieLine>
  gridpack::rtpr::RTPRDriver::getTieLines(
      gridpack::utility::Configuration::ChildCursors &tielines)
{
  // The ChildCursors argument is a vector of configuration
  // pointers. Each element in the vector is pointing at a seperate tieLine
  // block within the TieLines block in the input file.
  std::vector<gridpack::rtpr::TieLine> ret;
  int size = tielines.size();
  int i, idx;
  // Create string utilities object to help parse file
  gridpack::utility::StringUtils utils;
  // Loop over all child cursors
  for (idx = 0; idx < size; idx++) {
    std::string branchIDs;
    tielines[idx]->get("Branch",&branchIDs);
    std::vector<std::string> from_to = utils.blankTokenizer(branchIDs);
    gridpack::rtpr::TieLine tieline;
    tieline.from = atoi(from_to[0].c_str());
    tieline.to = atoi(from_to[1].c_str());
    std::string tag;
    tielines[idx]->get("Branch",&tag);
    std::string clean_tag = utils.clean2Char(tag);
    strncpy(tieline.tag,clean_tag.c_str(),2);
    tieline.tag[2] = '\0';
    ret.push_back(tieline);
  }
  return ret;
}

/**
 * Scale generation in a specified area
 * @param scale value to scale real power generation
 * @param area index of area
 */
void gridpack::rtpr::RTPRDriver::scaleAreaGeneration(double scale, int area)
{
  int nbus = p_pf_network->numBuses(); 
  int i;
  for (i=0; i<nbus; i++) {
    gridpack::powerflow::PFBus *bus = p_pf_network->getBus(i).get();
    if (bus->getArea() == area) {
      bus->scaleGeneratorRealPower(scale);
    }
  }
}

/**
 * Scale loads in a specified area
 * @param scale value to scale real power load
 * @param area index of area
 */
void gridpack::rtpr::RTPRDriver::scaleAreaLoads(double scale, int area)
{
  int nbus = p_pf_network->numBuses(); 
  int i;
  for (i=0; i<nbus; i++) {
    gridpack::powerflow::PFBus *bus = p_pf_network->getBus(i).get();
    if (bus->getArea() == area) {
      bus->scaleLoadRealPower(scale);
    }
  }
}

/**
 * Execute application. argc and argv are standard runtime parameters
 */
void gridpack::rtpr::RTPRDriver::execute(int argc, char** argv)
{
  int i, j;
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
  cursor = config->getCursor("Configuration.RealTimePathRating");
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
  int srcArea, dstArea;
  if (!cursor->get("sourceArea", &srcArea)) {
    srcArea = 1;
  }
  if (!cursor->get("destinationArea", &dstArea)) {
    dstArea = 1;
  }
  bool calcGenCntngcy;
  if (!cursor->get("calculateGeneratorContingencies",&tmp_bool)) {
    calcGenCntngcy = false;
  } else {
    util.toLower(tmp_bool);
    if (tmp_bool == "false") {
      calcGenCntngcy = false;
    } else {
      calcGenCntngcy = true;
    }
  }
  bool calcLineCntngcy;
  if (!cursor->get("calculateLineContingencies",&tmp_bool)) {
    calcLineCntngcy = false;
  } else {
    util.toLower(tmp_bool);
    if (tmp_bool == "false") {
      calcLineCntngcy = false;
    } else {
      calcLineCntngcy = true;
    }
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
  // Check for tie lines Set cursor so that it points to the
  // TieLines block in the input file
  cursor = config->getCursor(
      "ContingencyList.RealTimePathRating.tieLines");
  gridpack::utility::Configuration::ChildCursors
    tielines;
  if (cursor) cursor->children(tielines);
  std::vector<gridpack::rtpr::TieLine>  ties = getTieLines(tielines);
  int tsize = ties.size();
  std::vector<int> from_bus, to_bus;
  std::vector<std::string> tags;
  std::vector<bool> violations;
  from_bus.resize(tsize);
  to_bus.resize(tsize);
  tags.resize(tsize);
  for (i=0; i<tsize; i++) {
    from_bus[i] = ties[i].from;
    to_bus[i] = ties[i].to;
    tags[i] = ties[i].tag;
  }
  gridpack::parallel::Communicator task_comm = world.divide(grp_size);

  // Keep track of failed calculations
  std::vector<int> contingency_idx;
  std::vector<bool> contingency_success;
  // Keep track of which calculations completed successfully
  gridpack::parallel::GlobalVector<bool> ca_success(world);
  // Keep track of violation status for completed calculations
  // 1: no violations
  // 2: voltage violation
  // 3: line overload violation
  // 4: both voltage violation and line overload violation
  std::vector<int> contingency_violation;
  gridpack::parallel::GlobalVector<int> ca_violation(world);

  // Create powerflow applications on each task communicator
  boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network(new gridpack::powerflow::PFNetwork(task_comm));
  p_pf_network = pf_network;
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

  std::vector<gridpack::powerflow::Contingency> events;
  printf("p[%d] Got to 1\n",world.rank());
  if (!calcGenCntngcy && !calcLineCntngcy) {
    // Read in contingency file name
  printf("p[%d] Got to 2\n",world.rank());
    std::string contingencyfile;
    cursor = config->getCursor(
        "ContingencyList.RealTimePathRating");
    if (!cursor->get("contingencyList",&contingencyfile)) {
      contingencyfile = "contingencies.xml";
    }
  printf("p[%d] Got to 3\n",world.rank());
    if (world.rank() == 0) {
      printf("Contingency File: %s\n",contingencyfile.c_str());
    }
    // Open contingency file
    bool ok = config->open(contingencyfile,world);
  printf("p[%d] Got to 4\n",world.rank());

    // Get a list of contingencies. Set cursor so that it points to the
    // Contingencies block in the contingency file
    cursor = config->getCursor(
        "ContingencyList.RealTimePathRating.Contingencies");
    gridpack::utility::Configuration::ChildCursors contingencies;
  printf("p[%d] Got to 5\n",world.rank());
    if (cursor) cursor->children(contingencies);
    events = getContingencies(contingencies);
  } else {
  printf("p[%d] Got to 6\n",world.rank());
    std::vector<gridpack::powerflow::Contingency> genContingencies;
    if (calcGenCntngcy) {
      if (world.rank()==0) 
        printf("Calculating generator contingencies automatically\n");
      genContingencies = 
        createGeneratorContingencies(pf_network, srcArea);
    }
  printf("p[%d] Got to 7\n",world.rank());
    std::vector<gridpack::powerflow::Contingency> branchContingencies;
    if (calcLineCntngcy) {
      if (world.rank()==0) 
        printf("Calculating line contingencies automatically\n");
      branchContingencies =
        createBranchContingencies(pf_network, srcArea);
    }
    events = genContingencies;
    int nsize = branchContingencies.size();
    int ic;
    for (ic = 0; ic<nsize; ic++) {
      events.push_back(branchContingencies[ic]);
    }
  printf("p[%d] Got to 8\n",world.rank());
  }
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
      if (check_Qlim) pf_app.clearQlimViolations();
    } else {
#ifdef USE_SUCCESS
      contingency_success.push_back(false);
      contingency_violation.push_back(0);
#endif
      sprintf(sbuf,"\nDivergent for contingency %s\n",
          events[task_id].p_name.c_str());
      if (print_calcs) pf_app.print(sbuf);
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

  timer->stop(t_total);
  // If all processors executed at least one task, then print out timing
  // statistics (this printout does not work if some processors do not define
  // all timing variables)
  if (events.size()*grp_size >= world.size()) {
    timer->dump();
  }
}

