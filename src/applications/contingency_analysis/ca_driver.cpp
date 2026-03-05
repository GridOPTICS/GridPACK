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
 * @updated Yousu Chen
 * - N-1 auto-generation for branch and generator contingencies
 * - Automatic slack bus transfer and capacity check
 * - Q-limit support integration
 * @date  2026-01-31
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
#include "gridpack/utilities/results_exporter.hpp"
#include "ca_driver.hpp"

#include <sstream>

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
 * Auto-generate N-1 contingencies from the network
 * @param pf_app power flow application module with loaded network
 * @param gen_branches generate branch contingencies
 * @param gen_generators generate generator contingencies
 * @return vector of auto-generated contingencies
 */
std::vector<gridpack::powerflow::Contingency>
  gridpack::contingency_analysis::CADriver::generateN1Contingencies(
      gridpack::powerflow::PFAppModule &pf_app,
      bool gen_branches, bool gen_generators)
{
  std::vector<gridpack::powerflow::Contingency> ret;
  gridpack::utility::StringUtils utils;

  // Generate N-1 branch contingencies
  if (gen_branches) {
    std::vector<std::string> branch_data = pf_app.writeBranchString("flow_str");
    int branch_count = 0;

    for (size_t i=0; i<branch_data.size(); i++) {
      std::vector<std::string> tokens = utils.blankTokenizer(branch_data[i]);
      if (tokens.size()%8 != 0) {
        continue;  // Skip malformed data
      }

      int nline = tokens.size()/8;
      for (int j=0; j<nline; j++) {
        int from_bus = atoi(tokens[j*8].c_str());
        int to_bus = atoi(tokens[j*8+1].c_str());
        std::string ckt_id = tokens[j*8+2];

        // Create contingency for this branch
        gridpack::powerflow::Contingency contingency;
        char name_buf[64];
        sprintf(name_buf, "BR_%d_%d_%s", from_bus, to_bus,
                utils.clean2Char(ckt_id).c_str());
        contingency.p_name = name_buf;
        contingency.p_type = Branch;
        contingency.p_from.push_back(from_bus);
        contingency.p_to.push_back(to_bus);
        contingency.p_ckt.push_back(utils.clean2Char(ckt_id));
        contingency.p_saveLineStatus.push_back(true);

        ret.push_back(contingency);
        branch_count++;
      }
    }

    printf("Auto-generated %d N-1 branch contingencies\n", branch_count);
  }

  // Generate N-1 generator contingencies
  if (gen_generators) {
    std::vector<std::string> gen_data = pf_app.writeBusString("power");
    int gen_count = 0;

    for (size_t i=0; i<gen_data.size(); i++) {
      std::vector<std::string> tokens = utils.blankTokenizer(gen_data[i]);
      if (tokens.size()%4 != 0) {
        continue;  // Skip malformed data
      }

      int ngen = tokens.size()/4;
      for (int j=0; j<ngen; j++) {
        int bus_id = atoi(tokens[j*4].c_str());
        std::string gen_id = tokens[j*4+1];

        // Create contingency for this generator
        gridpack::powerflow::Contingency contingency;
        char name_buf[64];
        sprintf(name_buf, "GN_%d_%s", bus_id,
                utils.clean2Char(gen_id).c_str());
        contingency.p_name = name_buf;
        contingency.p_type = Generator;
        contingency.p_busid.push_back(bus_id);
        contingency.p_genid.push_back(utils.clean2Char(gen_id));
        contingency.p_saveGenStatus.push_back(true);

        ret.push_back(contingency);
        gen_count++;
      }
    }

    printf("Auto-generated %d N-1 generator contingencies\n", gen_count);
  }

  return ret;
}

/**
 * Check if a contingency is a duplicate of any in the existing list
 * @param contingency the contingency to check
 * @param existing_list vector of existing contingencies
 * @return true if duplicate found, false otherwise
 */
bool gridpack::contingency_analysis::CADriver::isDuplicateContingency(
    const gridpack::powerflow::Contingency &contingency,
    const std::vector<gridpack::powerflow::Contingency> &existing_list)
{
  for (size_t i = 0; i < existing_list.size(); i++) {
    const gridpack::powerflow::Contingency &existing = existing_list[i];

    // Check if same type
    if (existing.p_type != contingency.p_type) {
      continue;
    }

    if (contingency.p_type == Branch) {
      // For branch contingencies, check if same branches are tripped
      // A duplicate means all branches match (order doesn't matter)
      if (existing.p_from.size() != contingency.p_from.size()) {
        continue;
      }

      bool all_match = true;
      for (size_t j = 0; j < contingency.p_from.size(); j++) {
        bool found = false;
        for (size_t k = 0; k < existing.p_from.size(); k++) {
          if (contingency.p_from[j] == existing.p_from[k] &&
              contingency.p_to[j] == existing.p_to[k] &&
              contingency.p_ckt[j] == existing.p_ckt[k]) {
            found = true;
            break;
          }
        }
        if (!found) {
          all_match = false;
          break;
        }
      }

      if (all_match) {
        return true;  // Duplicate found
      }

    } else if (contingency.p_type == Generator) {
      // For generator contingencies, check if same generators are tripped
      if (existing.p_busid.size() != contingency.p_busid.size()) {
        continue;
      }

      bool all_match = true;
      for (size_t j = 0; j < contingency.p_busid.size(); j++) {
        bool found = false;
        for (size_t k = 0; k < existing.p_busid.size(); k++) {
          if (contingency.p_busid[j] == existing.p_busid[k] &&
              contingency.p_genid[j] == existing.p_genid[k]) {
            found = true;
            break;
          }
        }
        if (!found) {
          all_match = false;
          break;
        }
      }

      if (all_match) {
        return true;  // Duplicate found
      }
    }
  }

  return false;  // Not a duplicate
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
  // Check for Q limit violations (qlim: true=enabled, false=disabled)
  bool check_Qlim = cursor->get("qlim", true);
  // Output format: "json", "csv", or "text" (default)
  std::string outputFormat = "text";
  cursor->get("outputFormat", &outputFormat);
  std::string outputFile = "ca_results";
  cursor->get("outputFile", &outputFile);
  // Set static flag for PFBus class BEFORE network creation.
  // This controls how Q values are reported in output functions:
  // - When check_Qlim = false: output uses calculated Q from p_Qinj
  // - When check_Qlim = true: output uses p_qg (set by chkQlim())
  gridpack::powerflow::PFBus::setQlim(check_Qlim);
  gridpack::parallel::Communicator task_comm = world.divide(grp_size);

  // Keep track of failed calculations
#ifdef USE_SUCCESS
  std::vector<int> contingency_idx;
  std::vector<bool> contingency_success;
  gridpack::parallel::GlobalVector<bool> ca_success(world);
  std::vector<int> contingency_violation;
  gridpack::parallel::GlobalVector<int> ca_violation(world);
  std::vector<bool> contingency_isolated;
  gridpack::parallel::GlobalVector<bool> ca_isolated(world);
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

  // Collect base case results for export
  gridpack::utility::PowerFlowResults baseCaseResults;
  if (outputFormat != "text") {
    baseCaseResults = pf_app.collectResults();
  }

  // Check if auto-generation of N-1 contingencies is enabled
  // FullBranchN1: generate N-1 contingencies for all branches
  // FullGeneratorN1: generate N-1 contingencies for all generators
  cursor = config->getCursor("Configuration.Contingency_analysis");
  bool full_branch_n1 = false;
  bool full_generator_n1 = false;

  if (!cursor->get("FullBranchN1",&tmp_bool)) {
    full_branch_n1 = false;
  } else {
    util.toLower(tmp_bool);
    full_branch_n1 = (tmp_bool == "true");
  }

  if (!cursor->get("FullGeneratorN1",&tmp_bool)) {
    full_generator_n1 = false;
  } else {
    util.toLower(tmp_bool);
    full_generator_n1 = (tmp_bool == "true");
  }

  bool auto_generate_n1 = full_branch_n1 || full_generator_n1;

  std::vector<gridpack::powerflow::Contingency> events;
  int auto_generated_count = 0;
  int file_loaded_count = 0;
  int duplicates_skipped = 0;

  // Step 1: Auto-generate N-1 contingencies if requested
  if (auto_generate_n1) {
    if (world.rank() == 0) {
      printf("\n==================================================================\n");
      printf("Auto-generating N-1 contingencies from network\n");
      printf("  FullBranchN1: %s\n", full_branch_n1 ? "YES" : "NO");
      printf("  FullGeneratorN1: %s\n", full_generator_n1 ? "YES" : "NO");
      printf("==================================================================\n\n");
    }
    events = generateN1Contingencies(pf_app, full_branch_n1, full_generator_n1);
    auto_generated_count = events.size();
  }

  // Step 2: Load contingencies from file if specified
  // This allows combining auto-generated N-1 with custom N-2+ contingencies
  std::string contingencyfile;
  bool has_contingency_file = cursor->get("contingencyList",&contingencyfile);

  if (has_contingency_file || !auto_generate_n1) {
    // Set default filename if not specified
    if (!has_contingency_file) {
      contingencyfile = "contingencies.xml";
    }

    if (world.rank() == 0) {
      if (auto_generate_n1) {
        printf("Loading additional contingencies from file: %s\n", contingencyfile.c_str());
        printf("(Duplicates of auto-generated contingencies will be skipped)\n\n");
      } else {
        printf("Contingency List: %s\n", contingencyfile.c_str());
      }
    }

    // Open contingency file
    bool ok = config->open(contingencyfile,world);

    if (ok) {
      // Get a list of contingencies from file
      cursor = config->getCursor(
          "ContingencyList.Contingency_analysis.Contingencies");
      gridpack::utility::Configuration::ChildCursors contingencies;
      if (cursor) cursor->children(contingencies);
      std::vector<gridpack::powerflow::Contingency> file_contingencies =
          getContingencies(contingencies);

      // If auto-generation was used, check for duplicates before adding
      if (auto_generate_n1) {
        for (size_t i = 0; i < file_contingencies.size(); i++) {
          if (isDuplicateContingency(file_contingencies[i], events)) {
            duplicates_skipped++;
            if (world.rank() == 0) {
              printf("  Skipping duplicate: %s\n", file_contingencies[i].p_name.c_str());
            }
          } else {
            events.push_back(file_contingencies[i]);
            file_loaded_count++;
          }
        }
        if (world.rank() == 0 && file_loaded_count > 0) {
          printf("\nAdded %d unique contingencies from file\n", file_loaded_count);
          if (duplicates_skipped > 0) {
            printf("Skipped %d duplicates\n", duplicates_skipped);
          }
        }
      } else {
        // No auto-generation, just use file contingencies
        events = file_contingencies;
        file_loaded_count = events.size();
      }
    }
  }

  // Print summary
  if (world.rank() == 0) {
    printf("\n==================================================================\n");
    printf("Total contingencies to analyze: %d\n", (int)events.size());
    if (auto_generate_n1) {
      printf("  Auto-generated: %d\n", auto_generated_count);
      if (file_loaded_count > 0) {
        printf("  From file: %d\n", file_loaded_count);
      }
      if (duplicates_skipped > 0) {
        printf("  Duplicates skipped: %d\n", duplicates_skipped);
      }
    }
    printf("==================================================================\n\n");
  }

  // Print contingency details
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
  // Clear any Q limit warnings from base case before starting contingencies
  gridpack::powerflow::PFBus::clearQlimWarnings();

  // Local contingency results storage for JSON/CSV export
  std::vector<gridpack::utility::ContingencyResult> localContingencies;

  // Evaluate contingencies using the task manager
  int task_id;
  char sbuf[128];
  // nextTask returns the same task_id on all processors in task_comm. When the
  // calculation runs out of task, nextTask will return false.
  while (taskmgr.nextTask(task_comm, &task_id)) {
    printf("Executing task %d on process %d\n",task_id,world.rank());
    // Trim trailing spaces from contingency name for filename
    std::string fname = events[task_id].p_name;
    size_t end = fname.find_last_not_of(' ');
    if (end != std::string::npos) fname = fname.substr(0, end + 1);
    sprintf(sbuf,"%s.out",fname.c_str());
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
    // Sync ghost bus data after voltage reset to ensure branches connected to
    // ghost buses use the correct reset voltages in power flow calculation
    pf_network->updateBuses();
    // Set contingency
    bool contingencyFound = pf_app.setContingency(events[task_id]);
    if (!contingencyFound) {
      printf("WARNING: Contingency '%s' - elements not found or no valid slack bus\n",
             events[task_id].p_name.c_str());
    }
    // Check for islanding before attempting to solve
    // Note: lone bus isolation is handled separately as a warning, not a failure
    int islandCount = pf_app.getIslandCount();
    bool hasLoneBus = pf_app.hasLoneBus();
    bool islandDetected = (islandCount > 1);
    // Solve power flow equations for this system
#ifdef USE_SUCCESS
    contingency_idx.push_back(task_id);
#endif
    // Skip power flow if contingency setup failed (no valid slack) or islanding detected
    bool slackCapacityOk = true;  // Will be checked after solve
    if (contingencyFound && !islandDetected && pf_app.solve()) {
      if (check_Qlim && !pf_app.checkQlimViolations()) {
        pf_app.solve();
      }
      // Write PV->PQ conversion warnings to output file
      if (print_calcs && check_Qlim) {
        std::vector<std::string>& warnings = gridpack::powerflow::PFBus::getQlimWarnings();
        for (size_t w = 0; w < warnings.size(); w++) {
          pf_app.print(warnings[w].c_str());
        }
      }
      // Check if slack bus generator exceeds capacity
      slackCapacityOk = pf_app.checkSlackCapacity();
      if (!slackCapacityOk) {
        // Slack generator exceeds Pmax - insufficient generation capacity
        // This is treated as a failure, similar to divergence
#ifdef USE_SUCCESS
        contingency_success.push_back(false);
        contingency_violation.push_back(0);
        contingency_isolated.push_back(false);
#endif
        if (outputFormat != "text") {
          gridpack::utility::ContingencyResult ctResult;
          ctResult.name = events[task_id].p_name;
          ctResult.type = (events[task_id].p_type == Branch) ? "branch" : "generator";
          ctResult.hasVoltageViolation = false;
          ctResult.hasBranchViolation = false;
          ctResult.solution.convergence.converged = false;
          localContingencies.push_back(ctResult);
        }
        sprintf(sbuf,"\nInsufficient generation capacity for contingency %s\n",
            events[task_id].p_name.c_str());
        if (print_calcs) pf_app.print(sbuf);
      } else {
        // Power flow solved and slack within capacity
#ifdef USE_SUCCESS
        contingency_success.push_back(true);
        contingency_isolated.push_back(hasLoneBus);
#endif
        // If power flow solution is successful, write out voltages and currents
        if (print_calcs) pf_app.write();
        // Check for violations
        bool ok1 = pf_app.checkVoltageViolations();
        bool ok2 = pf_app.checkLineOverloadViolations();
        bool ok = ok1 && ok2;
        // Collect results for JSON/CSV export
        if (outputFormat != "text") {
          gridpack::utility::ContingencyResult ctResult;
          ctResult.name = events[task_id].p_name;
          ctResult.type = (events[task_id].p_type == Branch) ? "branch" : "generator";
          ctResult.hasVoltageViolation = !ok1;
          ctResult.hasBranchViolation = !ok2;
          ctResult.solution = pf_app.collectResults();
          localContingencies.push_back(ctResult);
        }
      // Include results of violation checks in output
      if (ok) {
        sprintf(sbuf,"\nNo violation for contingency %s\n",
            events[task_id].p_name.c_str());
#ifdef USE_SUCCESS
        contingency_violation.push_back(1);
#endif
      }
      // Report bus voltage violations
      if (!ok1) {
        sprintf(sbuf,"\nBus Violation for contingency %s\n",
            events[task_id].p_name.c_str());
      } else if (!ok) {
        sprintf(sbuf,"\nNo Bus Violation for contingency %s\n",
            events[task_id].p_name.c_str());
      }
      if (print_calcs) pf_app.print(sbuf);
      if (print_calcs) pf_app.writeCABus();
      // Report branch overload violations
      if (!ok2) {
        sprintf(sbuf,"\nBranch Violation for contingency %s\n",
            events[task_id].p_name.c_str());
      } else if (!ok) {
        sprintf(sbuf,"\nNo Branch Violation for contingency %s\n",
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
        // Note: clearQlimViolations() moved after unSetContingency() below
      }  // end slackCapacityOk block
    } else {
#ifdef USE_SUCCESS
      contingency_success.push_back(false);
      contingency_violation.push_back(0);
      contingency_isolated.push_back(false);
#endif
      if (outputFormat != "text") {
        gridpack::utility::ContingencyResult ctResult;
        ctResult.name = events[task_id].p_name;
        ctResult.type = (events[task_id].p_type == Branch) ? "branch" : "generator";
        ctResult.hasVoltageViolation = false;
        ctResult.hasBranchViolation = false;
        ctResult.solution.convergence.converged = false;
        localContingencies.push_back(ctResult);
      }
      if (islandDetected) {
        sprintf(sbuf,"\nIslanding detected for contingency %s (%d islands)\n",
            events[task_id].p_name.c_str(), islandCount);
      } else if (!contingencyFound) {
        sprintf(sbuf,"\nNo valid slack bus for contingency %s\n",
            events[task_id].p_name.c_str());
      } else {
        sprintf(sbuf,"\nDivergent for contingency %s\n",
            events[task_id].p_name.c_str());
      }
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
    // Clear Q limit violations AFTER unSetContingency so generators are restored first.
    // This ensures clearQlim() sees the correct generator status when deciding
    // whether to restore p_isPV (PV bus status).
    if (check_Qlim) pf_app.clearQlimViolations();
    // Clear Q limit warnings for next contingency
    gridpack::powerflow::PFBus::clearQlimWarnings();
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
    ca_isolated.addElements(contingency_idx, contingency_isolated);
  }
  ca_success.upload();
  ca_violation.upload();
  ca_isolated.upload();
  // All processes call getData to ensure GA progress (NGA_Gather requires
  // remote process participation for one-sided communication).
  contingency_idx.clear();
  contingency_success.clear();
  contingency_violation.clear();
  contingency_isolated.clear();
  for (i=0; i<ntasks; i++) contingency_idx.push_back(i);
  ca_success.getData(contingency_idx, contingency_success);
  contingency_success.clear();
  ca_violation.getData(contingency_idx, contingency_violation);
  ca_isolated.getData(contingency_idx, contingency_isolated);
  // Write out stats on successful calculations
  if (world.rank() == 0) {
    std::ofstream fout;
    fout.open("success.txt");
    for (i=0; i<ntasks; i++) {
      if (contingency_success[i]) {
        fout << "contingency: " << i+1 << " success: true";
        if (contingency_violation[i] == 1) {
          fout << " violation: none";
        } else if (contingency_violation[i] == 2) {
          fout << " violation: bus";
        } else if (contingency_violation[i] == 3) {
          fout << " violation: branch";
        } else if (contingency_violation[i] == 4) {
          fout << " violation: bus and branch";
        }
        if (contingency_isolated[i]) {
          fout << " warning: isolated";
        }
        fout << std::endl;
      } else {
        fout << "contingency: " << i+1 << " success: false" << std::endl;
      }
    }
    fout.close();
  }
#endif

  // Export CA results to JSON or CSV
  if (outputFormat == "csv") {
    // Rank 0 writes base case (creates files with headers)
    if (world.rank() == 0) {
      gridpack::utility::ResultsExporter::writePFCSV(outputFile,
          baseCaseResults, "base_case");
    }
    world.barrier();
    // Each process writes its contingencies in order (files opened in append mode)
    for (int p = 0; p < world.size(); p++) {
      if (world.rank() == p) {
        for (size_t ci = 0; ci < localContingencies.size(); ci++) {
          gridpack::utility::ResultsExporter::writePFCSV(
              outputFile, localContingencies[ci].solution,
              localContingencies[ci].name);
        }
      }
      world.barrier();
    }
    // Rank 0 writes summary CSV using gathered data
#ifdef USE_SUCCESS
    if (world.rank() == 0) {
      std::string summaryFile = outputFile + "_summary.csv";
      std::ofstream sout(summaryFile.c_str());
      sout << "contingency,type,converged,has_voltage_violation,has_branch_violation\n";
      for (int ci = 0; ci < ntasks; ci++) {
        bool converged = (contingency_violation[ci] > 0);
        sout << events[ci].p_name << ","
             << (events[ci].p_type == Branch ? "branch" : "generator") << ","
             << (converged ? "true" : "false") << ","
             << ((contingency_violation[ci] == 2 || contingency_violation[ci] == 4)
                 ? "true" : "false") << ","
             << ((contingency_violation[ci] == 3 || contingency_violation[ci] == 4)
                 ? "true" : "false") << "\n";
      }
      sout.close();
    }
#endif
  }

  if (outputFormat == "json") {
    // Each process serializes its contingency results as JSON text
    std::ostringstream localJsonStream;
    for (size_t ci = 0; ci < localContingencies.size(); ci++) {
      gridpack::utility::ResultsExporter::writeContingencyResultJSON(
          localJsonStream, localContingencies[ci]);
      if (ci + 1 < localContingencies.size()) {
        localJsonStream << ",\n";
      }
    }
    std::string localJson = localJsonStream.str();

    // Gather JSON fragment sizes
    int localLen = static_cast<int>(localJson.size());
    std::vector<int> allLens(world.size());
    MPI_Allgather(&localLen, 1, MPI_INT, &allLens[0], 1, MPI_INT,
                  static_cast<MPI_Comm>(world));

    // Compute displacements
    std::vector<int> displs(world.size());
    int totalLen = 0;
    for (int p = 0; p < world.size(); p++) {
      displs[p] = totalLen;
      totalLen += allLens[p];
    }

    // Gather all JSON fragments on rank 0
    std::vector<char> allJsonBuf;
    if (world.rank() == 0) allJsonBuf.resize(totalLen);
    MPI_Gatherv(localJson.c_str(), localLen, MPI_CHAR,
                world.rank() == 0 ? &allJsonBuf[0] : NULL,
                &allLens[0], &displs[0], MPI_CHAR,
                0, static_cast<MPI_Comm>(world));

    // Rank 0 writes the final JSON file
    if (world.rank() == 0) {
      std::string jsonFile = outputFile + ".json";
      std::ofstream jout(jsonFile.c_str());
      // Write base case and open contingencies array
      gridpack::utility::ResultsExporter::writeCAJSONHeader(jout,
          baseCaseResults);
      // Write contingencies from gathered fragments
      bool firstFragment = true;
      for (int p = 0; p < world.size(); p++) {
        if (allLens[p] > 0) {
          if (!firstFragment) jout << ",\n";
          jout.write(&allJsonBuf[displs[p]], allLens[p]);
          firstFragment = false;
        }
      }
      jout << "\n";
      gridpack::utility::ResultsExporter::writeCAJSONFooter(jout);
      jout.close();
    }
  }

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
  if (events.size()*grp_size >= world.size()) {
    timer->dump();
  }
}

