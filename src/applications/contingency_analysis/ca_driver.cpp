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
 * @updated Yousu Chen
 * - csv_flat / csv_delta per-(contingency,branch) outputs
 * - monitorBranchesFile / monitorAreas / monitorKvMin/Max filters (all formats)
 * @date  2026-06-21
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

#include <boost/scoped_ptr.hpp>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>

// Statistical-summary output (vmag.txt, pflow.txt, etc.) used to be controlled
// by a USE_STATBLOCK build-time macro; it is now a runtime XML option,
// `Configuration.Contingency_analysis.writeStats`, defaulting to true to
// preserve existing behavior.
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
      if (!contingencies[idx]->get("CKT",&names)) {
        contingencies[idx]->get("contingencyLineNames",&names);
      }
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
      if (!contingencies[idx]->get("GenID",&gens)) {
        contingencies[idx]->get("contingencyGenerators",&gens);
      }
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

    if (gridpack::parallel::Communicator().rank() == 0) {
      printf("Auto-generated %d N-1 branch contingencies\n", branch_count);
    }
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

    if (gridpack::parallel::Communicator().rank() == 0) {
      printf("Auto-generated %d N-1 generator contingencies\n", gen_count);
    }
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
  // Statistical-summary output (vmag.txt, pflow.txt, etc. via StatBlock).
  // Default true to preserve existing behavior; set false to skip the
  // per-case StatBlock work and the 13 post-loop global writes.
  bool write_stats = true;
  if (cursor->get("writeStats",&tmp_bool)) {
    util.toLower(tmp_bool);
    write_stats = (tmp_bool != "false");
  }
  // groupSize is forced to 1: an outaged branch may straddle a multi-rank
  // partition. Scale by adding ranks, not by widening a group.
  if (!cursor->get("groupSize",&grp_size)) {
    grp_size = 1;
  }
  if (grp_size != 1) {
    if (world.rank() == 0) {
      printf("WARNING: groupSize=%d is not supported (a contingency branch "
             "may span the partition boundary); using groupSize=1\n",
             grp_size);
    }
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
  double qlim_deadband = cursor->get("qlimDeadband", 0.1);
  // Output format: "text" (default), "json", "csv", "csv_flat", "csv_delta".
  std::string outputFormat = "text";
  cursor->get("outputFormat", &outputFormat);
  if (outputFormat != "text" && outputFormat != "json" &&
      outputFormat != "csv"  && outputFormat != "csv_flat" &&
      outputFormat != "csv_delta") {
    if (world.rank() == 0) {
      printf("ERROR: unrecognized outputFormat='%s'. "
             "Must be one of: text, json, csv, csv_flat, csv_delta. Aborting.\n",
             outputFormat.c_str());
    }
    world.barrier();
    MPI_Abort(static_cast<MPI_Comm>(world), 1);
  }
  std::string outputFile = "ca_results";
  cursor->get("outputFile", &outputFile);
  // Optional CSV allowlist (from_bus,to_bus,ckt). Empty -> emit all.
  std::string monitorBranchesFile;
  cursor->get("monitorBranchesFile", &monitorBranchesFile);
  // Optional area/kV gates. Empty/zero/missing -> no restriction on that
  // dimension. Filters AND together with monitorBranchesFile.
  std::string monitorAreasStr;
  cursor->get("monitorAreas", &monitorAreasStr);
  double monitorKvMin = 0.0;
  cursor->get("monitorKvMin", &monitorKvMin);
  double monitorKvMax = 0.0;
  cursor->get("monitorKvMax", &monitorKvMax);
  std::set<int> monitorAreas;
  {
    std::vector<std::string> tok = util.blankTokenizer(monitorAreasStr);
    for (size_t i = 0; i < tok.size(); i++) {
      if (!tok[i].empty()) monitorAreas.insert(atoi(tok[i].c_str()));
    }
  }
  // Any monitor filter configured; gates every output format.
  bool haveMonitorFilter = !monitorAreas.empty() ||
                           monitorKvMin > 0.0 || monitorKvMax > 0.0 ||
                           !monitorBranchesFile.empty();
  // Loading% denominator for all CA outputs (.out, _violations.csv, JSON,
  // csv_flat, csv_delta). A|B|C, default A to match PW/PSSE convention.
  // A->B->C fallback if the requested tier is zero/missing.
  std::string contingencyRating = "A";
  cursor->get("contingencyRating", &contingencyRating);
  util.toUpper(contingencyRating);
  if (contingencyRating != "A" && contingencyRating != "B" &&
      contingencyRating != "C") {
    if (world.rank() == 0) {
      printf("WARNING: contingencyRating='%s' not A/B/C; defaulting to A\n",
             contingencyRating.c_str());
    }
    contingencyRating = "A";
  }
  // Severity threshold for violation reporting; loading% > threshold*100
  // is flagged. Default 1.0 (100% of rate).
  double violationSeverityThreshold = 1.0;
  cursor->get("violationSeverityThreshold", &violationSeverityThreshold);
  if (violationSeverityThreshold <= 0.0) violationSeverityThreshold = 1.0;
  // Cap on top_severe_contingencies and roster arrays in _summary.json.
  int topN = 10;
  cursor->get("topN", &topN);
  if (topN < 1) topN = 1;
  if (topN > 10000) topN = 10000;
  // Weights for composite_pi = piBranchWeight*branch_pi + piVoltageWeight*voltage_pi.
  double piBranchWeight = 1.0, piVoltageWeight = 1.0;
  cursor->get("piBranchWeight",  &piBranchWeight);
  cursor->get("piVoltageWeight", &piVoltageWeight);
  if (piBranchWeight  < 0.0) piBranchWeight  = 0.0;
  if (piVoltageWeight < 0.0) piVoltageWeight = 0.0;

  // Set static flag for PFBus class BEFORE network creation.
  // This controls how Q values are reported in output functions:
  // - When check_Qlim = false: output uses calculated Q from p_Qinj
  // - When check_Qlim = true: output uses p_qg (set by chkQlim())
  gridpack::powerflow::PFBus::setQlim(check_Qlim);
  gridpack::powerflow::PFBus::setQlimDeadband(qlim_deadband);
  gridpack::parallel::Communicator task_comm = world.divide(grp_size);

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

  // Build (number -> name) lookup tables for area, zone, owner.
  // Keyed on the PSS/E-assigned number (not contiguous), used when
  // emitting per-(branch,contingency) CSV rows so each row carries
  // human-readable area/zone/owner names alongside the numbers.
  std::map<int, std::string> area_name_by_num;
  std::map<int, std::string> zone_name_by_num;
  std::map<int, std::string> owner_name_by_num;
  {
    boost::shared_ptr<gridpack::component::DataCollection> netdata =
      pf_network->getNetworkData();
    int aT = 0, zT = 0, oT = 0;
    netdata->getValue(AREA_TOTAL,  &aT);
    netdata->getValue(ZONE_TOTAL,  &zT);
    netdata->getValue(OWNER_TOTAL, &oT);
    for (int i = 0; i < aT; i++) {
      int n = 0; std::string s;
      netdata->getValue(AREAINTG_NUMBER, &n, i);
      netdata->getValue(AREAINTG_NAME,   &s, i);
      area_name_by_num[n] = s;
    }
    for (int i = 0; i < zT; i++) {
      int n = 0; std::string s;
      netdata->getValue(ZONE_NUMBER, &n, i);
      netdata->getValue(ZONE_NAME,   &s, i);
      zone_name_by_num[n] = s;
    }
    for (int i = 0; i < oT; i++) {
      int n = 0; std::string s;
      netdata->getValue(OWNER_NUMBER, &n, i);
      netdata->getValue(OWNER_NAME,   &s, i);
      owner_name_by_num[n] = s;
    }
  }

  // Per-rank bus metadata + wide/long branch-row outputs for csv_flat
  // (long-form one row per branch per case) and csv_delta (wide-form one
  // row per branch per case joining base and cont state). Both share the
  // bus_meta load, the buses sidecar, and the per-rank .part-file gather.
  bool wantBusSidecar = (outputFormat == "csv_flat" ||
                         outputFormat == "csv_delta");
  struct BusMeta {
    std::string name;
    double basekv;
    int    area, zone, owner;
  };
  std::map<int, BusMeta> bus_meta;
  if (wantBusSidecar) {
    int nBus = pf_network->numBuses();
    for (int i = 0; i < nBus; i++) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>(pf_network->getBus(i).get());
      if (!bus) continue;
      int orig = pf_network->getOriginalBusIndex(i);
      BusMeta m;
      m.name   = bus->getBusName();
      m.basekv = bus->getBaseKV();
      m.area   = bus->getArea();
      m.zone   = bus->getZone();
      m.owner  = bus->getOwner();
      bus_meta[orig] = m;
    }
  }
  // (area, base kV) per bus, all-gathered over task_comm so filter lookups
  // work on gathered rows.
  struct BusAreaKv { int area; double basekv; };
  std::map<int, BusAreaKv> bus_ak;
  if (wantBusSidecar || haveMonitorFilter) {
    std::vector<int> lid, larea;
    std::vector<double> lkv;
    int nBus = pf_network->numBuses();
    for (int i = 0; i < nBus; i++) {
      if (!pf_network->getActiveBus(i)) continue;
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>(pf_network->getBus(i).get());
      if (!bus) continue;
      lid.push_back(pf_network->getOriginalBusIndex(i));
      larea.push_back(bus->getArea());
      lkv.push_back(bus->getBaseKV());
    }
    MPI_Comm tc = static_cast<MPI_Comm>(task_comm);
    int tsize = task_comm.size();
    int nloc = static_cast<int>(lid.size());
    std::vector<int> counts(tsize, 0), displs(tsize, 0);
    MPI_Allgather(&nloc, 1, MPI_INT, &counts[0], 1, MPI_INT, tc);
    int tot = 0;
    for (int p = 0; p < tsize; p++) { displs[p] = tot; tot += counts[p]; }
    std::vector<int> gid(tot > 0 ? tot : 1), garea(tot > 0 ? tot : 1);
    std::vector<double> gkv(tot > 0 ? tot : 1);
    int *lid_p   = nloc > 0 ? &lid[0]   : NULL;
    int *larea_p = nloc > 0 ? &larea[0] : NULL;
    double *lkv_p = nloc > 0 ? &lkv[0]  : NULL;
    MPI_Allgatherv(lid_p,   nloc, MPI_INT,    &gid[0],   &counts[0], &displs[0], MPI_INT,    tc);
    MPI_Allgatherv(larea_p, nloc, MPI_INT,    &garea[0], &counts[0], &displs[0], MPI_INT,    tc);
    MPI_Allgatherv(lkv_p,   nloc, MPI_DOUBLE, &gkv[0],   &counts[0], &displs[0], MPI_DOUBLE, tc);
    for (int i = 0; i < tot; i++) {
      BusAreaKv a;
      a.area = garea[i];
      a.basekv = gkv[i];
      bus_ak[gid[i]] = a;
    }
  }

  // Strip surrounding single quotes (PSS/E style) and outer whitespace.
  auto trim_quoted = [](const std::string &in) -> std::string {
    std::string s = in;
    size_t a = s.find_first_not_of(" \t");
    size_t b = s.find_last_not_of(" \t");
    if (a == std::string::npos) return std::string();
    s = s.substr(a, b - a + 1);
    if (s.size() >= 2 && s.front() == '\'' && s.back() == '\'') {
      s = s.substr(1, s.size() - 2);
    }
    a = s.find_first_not_of(" \t");
    b = s.find_last_not_of(" \t");
    return (a == std::string::npos) ? std::string() : s.substr(a, b - a + 1);
  };
  auto lookup_name = [&](const std::map<int, std::string> &m, int n) -> std::string {
    std::map<int, std::string>::const_iterator it = m.find(n);
    if (it == m.end()) return std::string();
    return trim_quoted(it->second);
  };

  // Per-rank bus metadata sidecar (deduped by world rank 0 after the loop).
  if (wantBusSidecar) {
    std::ostringstream oss;
    oss << outputFile << "_buses." << world.rank() << ".part";
    std::ofstream fbus(oss.str().c_str(),
                       std::ios::out | std::ios::trunc | std::ios::binary);
    fbus << std::fixed;
    for (std::map<int, BusMeta>::const_iterator it = bus_meta.begin();
         it != bus_meta.end(); ++it) {
      const BusMeta &m = it->second;
      fbus << it->first << ","
           << trim_quoted(m.name) << ","
           << std::setprecision(2) << m.basekv << ","
           << m.area << "," << m.zone << "," << m.owner << ","
           << lookup_name(area_name_by_num,  m.area)  << ","
           << lookup_name(zone_name_by_num,  m.zone)  << ","
           << lookup_name(owner_name_by_num, m.owner)
           << "\n";
    }
    fbus.close();
  }

  // Per-rank streaming output file. Opened on first row written so non-csv_flat
  // runs and ranks that produce no rows leave nothing behind.
  std::string flatPartPath;
  if (outputFormat == "csv_flat") {
    std::ostringstream oss;
    oss << outputFile << "_flat." << world.rank() << ".part";
    flatPartPath = oss.str();
  }
  std::ofstream flatPart;
  size_t flatRowCount = 0;

  // Per-rank _violations.csv stream. Populated by every output format that
  // knows about violations (json/csv via populateViolations; csv_flat/csv_delta
  // inline where loading_pct is already computed). Emit-branch and emit-voltage
  // helpers below open the file lazily.
  std::string violPartPath;
  {
    std::ostringstream oss;
    oss << outputFile << "_violations." << world.rank() << ".part";
    violPartPath = oss.str();
  }
  std::ofstream violPart;
  size_t violRowCount = 0;
  // Running summary state; captureFlatRows / captureDeltaRows / populateViolations
  // all funnel through the emit helpers, so this stays consistent regardless of
  // outputFormat.
  struct WorstBranchState {
    double loading_pct = 0.0;
    int from = 0, to = 0;
    std::string ckt;
    std::string ct_name;
  };
  struct WorstVoltageState {
    double v_pu = 1.0;
    double dev_pu = 0.0;   // signed
    int bus_id = 0;
    std::string ct_name;
  };
  WorstBranchState worstBr;
  WorstVoltageState worstVLo;
  worstVLo.v_pu = 1e9;   // start high so any real low value beats it
  WorstVoltageState worstVHi;
  worstVHi.v_pu = -1e9;
  std::set<std::string> ctsWithBranchViol;
  std::set<std::string> ctsWithVoltageViol;
  // Rank-local per-ct composite indices; reduced to world 0 in summary.
  std::map<std::string, double> ctPi;             // sum(mva/rate)^2, all monitored branches
  std::map<std::string, double> ctVpi;            // sum((v-1)/dv)^2, all energized buses
  std::map<std::string, double> ctVdev;           // sum(v-limit)^2, violated buses only
  std::map<std::string, double> ctWorstLoading;   // max loading_pct among violated branches
  std::map<std::string, double> ctWorstVdev;      // max |v - limit| among violated buses
  std::map<std::string, double> ctWorstVpu;       // v_pu that produced ctWorstVdev
  auto accumBranchPi = [&](const std::string &ct_name,
                           double mva, double rate) {
    if (rate <= 0.0) return;
    double r = mva / rate;
    ctPi[ct_name] += r * r;
  };
  // Textbook voltage PI: normalize by the direction-aware half-band so a bus
  // at Vmin or Vmax contributes 1.0 even when limits are asymmetric.
  // Skip isolated (v<=0) and NaN/inf buses so numerical failures on a single
  // bus don't poison the accumulator.
  auto accumVoltagePi = [&](const std::string &ct_name, double v_pu) {
    if (v_pu <= 0.0 || !std::isfinite(v_pu)) return;
    double denom = (v_pu < 1.0) ? (1.0 - Vmin) : (Vmax - 1.0);
    if (denom <= 0.0) return;
    double r = (v_pu - 1.0) / denom;
    ctVpi[ct_name] += r * r;
  };
  // Legacy deviation metric: (v-limit)^2 accrued only for violated buses.
  auto accumVoltageDev = [&](const std::string &ct_name, double dev_pu) {
    ctVdev[ct_name] += dev_pu * dev_pu;
  };
  auto openViolPart = [&]() {
    if (violPart.is_open()) return;
    violPart.open(violPartPath.c_str(), std::ios::out | std::ios::trunc);
    violPart << std::fixed;
  };
  // Row schema (same 11 columns for branch and voltage rows; unused fields blank):
  //   event_idx,contingency,type,element,mva_or_vpu,rate_or_limit,loading_percent,
  //   base_mva,delta,severity
  // For branch: element = "from-to-ckt", mva_or_vpu = MVA,      rate_or_limit = rate,
  //             loading_percent = %,     base_mva/delta populated,        severity.
  // For voltage: element = "bus_id",     mva_or_vpu = v_pu,     rate_or_limit = "low/high",
  //             loading_percent = "",    base_mva = "",  delta = signed dev pu.
  auto emitBranchViolation = [&](int event_idx, const std::string &ct_name,
                                 int from, int to, const std::string &ckt,
                                 double mva, double rate,
                                 double loading_pct, double base_mva) {
    openViolPart();
    const char *sev = (loading_pct >= 105.0) ? "critical" : "warning";
    violPart << event_idx << "," << ct_name << ",branch,"
             << from << "-" << to << "-" << ckt << ","
             << std::setprecision(4) << mva << ","
             << std::setprecision(4) << rate << ","
             << std::setprecision(2) << loading_pct << ","
             << std::setprecision(4) << base_mva << ","
             << std::setprecision(4) << (mva - base_mva) << ","
             << sev << "\n";
    violRowCount++;
    if (loading_pct > worstBr.loading_pct) {
      worstBr.loading_pct = loading_pct;
      worstBr.from = from; worstBr.to = to;
      worstBr.ckt = ckt;
      worstBr.ct_name = ct_name;
    }
    ctsWithBranchViol.insert(ct_name);
    double &pc = ctWorstLoading[ct_name];
    if (loading_pct > pc) pc = loading_pct;
  };
  auto emitVoltageViolation = [&](int event_idx, const std::string &ct_name,
                                  int bus_id, double v_pu,
                                  double lo, double hi) {
    openViolPart();
    double dev = (v_pu < lo) ? (v_pu - lo)
               : (v_pu > hi) ? (v_pu - hi)
               : 0.0;
    if (dev == 0.0) return;
    const char *sev = (std::abs(dev) >= 0.05) ? "critical" : "warning";
    // rate_or_limit column holds "low_pu:high_pu" for voltage rows.
    std::ostringstream limits;
    limits << std::setprecision(4) << std::fixed << lo << ":" << hi;
    violPart << event_idx << "," << ct_name << ",voltage,"
             << bus_id << ","
             << std::setprecision(6) << v_pu << ","
             << limits.str() << ","
             << ","                                 // loading_percent blank
             << ","                                 // base_mva blank
             << std::setprecision(6) << dev << ","
             << sev << "\n";
    violRowCount++;
    if (dev < 0.0 && v_pu < worstVLo.v_pu) {
      worstVLo.v_pu = v_pu;
      worstVLo.dev_pu = dev;
      worstVLo.bus_id = bus_id;
      worstVLo.ct_name = ct_name;
    }
    if (dev > 0.0 && v_pu > worstVHi.v_pu) {
      worstVHi.v_pu = v_pu;
      worstVHi.dev_pu = dev;
      worstVHi.bus_id = bus_id;
      worstVHi.ct_name = ct_name;
    }
    ctsWithVoltageViol.insert(ct_name);
    accumVoltageDev(ct_name, dev);
    double absDev = std::abs(dev);
    double &d = ctWorstVdev[ct_name];
    if (absDev > d) { d = absDev; ctWorstVpu[ct_name] = v_pu; }
  };

  // (from, to, ckt) key shared by the monitor allowlist and base_cache.
  struct BranchKey {
    int from, to;
    std::string ckt;
    bool operator<(const BranchKey &o) const {
      if (from != o.from) return from < o.from;
      if (to   != o.to  ) return to   < o.to;
      return ckt < o.ckt;
    }
  };

  // Per-branch rate-A/B/C from the parsed network data. Keyed by (from,to,ckt)
  // so the csv_flat / csv_delta emit paths can pick the configured rating.
  // Built once before the contingency loop. Each rank only sees its own
  // active+ghost branches; that's fine -- the emit path is also rank-local.
  struct BranchRates {
    double rate_a, rate_b, rate_c;
  };
  std::map<BranchKey, BranchRates> branch_rates;
  if (outputFormat == "csv_flat" || outputFormat == "csv_delta") {
    int nBranch = pf_network->numBranches();
    for (int i = 0; i < nBranch; i++) {
      boost::shared_ptr<gridpack::component::DataCollection> bd =
        pf_network->getBranchData(i);
      if (!bd) continue;
      int from = 0, to = 0, nelems = 0;
      bd->getValue(BRANCH_FROMBUS, &from);
      bd->getValue(BRANCH_TOBUS,   &to);
      if (!bd->getValue(BRANCH_NUM_ELEMENTS, &nelems)) continue;
      for (int k = 0; k < nelems; k++) {
        std::string ckt;
        if (!bd->getValue(BRANCH_CKT, &ckt, k)) continue;
        // Trim leading/trailing whitespace and PSS/E surrounding quotes.
        size_t a = ckt.find_first_not_of(" \t");
        size_t b = ckt.find_last_not_of(" \t");
        ckt = (a == std::string::npos) ? std::string()
                                       : ckt.substr(a, b - a + 1);
        if (ckt.size() >= 2 && ckt.front() == '\'' && ckt.back() == '\'') {
          ckt = ckt.substr(1, ckt.size() - 2);
          a = ckt.find_first_not_of(" \t");
          b = ckt.find_last_not_of(" \t");
          ckt = (a == std::string::npos) ? std::string()
                                         : ckt.substr(a, b - a + 1);
        }
        BranchRates r;
        r.rate_a = 0.0; r.rate_b = 0.0; r.rate_c = 0.0;
        bd->getValue(BRANCH_RATING_A, &r.rate_a, k);
        bd->getValue(BRANCH_RATING_B, &r.rate_b, k);
        bd->getValue(BRANCH_RATING_C, &r.rate_c, k);
        BranchKey key;
        key.from = from; key.to = to; key.ckt = ckt;
        branch_rates[key] = r;
      }
    }
  }
  // Base case always uses rate-A (PSS/E "normal" rating). Contingency rows
  // use whichever the user picked, with A->B->C fallback if zero/missing.
  auto pickContRate = [&](const BranchRates &r) -> double {
    if (contingencyRating == "A") {
      return r.rate_a;
    }
    if (contingencyRating == "B") {
      return (r.rate_b > 0.0) ? r.rate_b : r.rate_a;
    }
    if (r.rate_c > 0.0) return r.rate_c;
    if (r.rate_b > 0.0) return r.rate_b;
    return r.rate_a;
  };

  // Monitor allowlist parsed from monitorBranchesFile. Empty -> emit all.
  std::set<BranchKey> monitorSet;
  if (!monitorBranchesFile.empty()) {
    std::ifstream fin(monitorBranchesFile.c_str());
    if (!fin.is_open()) {
      if (world.rank() == 0) {
        printf("WARNING: monitorBranchesFile '%s' not found; emitting all branches\n",
               monitorBranchesFile.c_str());
      }
    } else {
      std::string line;
      size_t lineNo = 0;
      while (std::getline(fin, line)) {
        lineNo++;
        // Strip trailing CR (Windows line endings).
        while (!line.empty() && (line[line.size()-1] == '\r' ||
                                 line[line.size()-1] == '\n')) {
          line.resize(line.size()-1);
        }
        // Skip blank lines and comments.
        size_t firstNon = line.find_first_not_of(" \t");
        if (firstNon == std::string::npos) continue;
        if (line[firstNon] == '#') continue;
        // Tokenize on commas.
        std::vector<std::string> tok;
        size_t pos = 0;
        while (pos <= line.size()) {
          size_t comma = line.find(',', pos);
          std::string t = (comma == std::string::npos)
                          ? line.substr(pos)
                          : line.substr(pos, comma - pos);
          size_t a = t.find_first_not_of(" \t");
          size_t b = t.find_last_not_of(" \t");
          tok.push_back((a == std::string::npos) ? std::string()
                                                 : t.substr(a, b - a + 1));
          if (comma == std::string::npos) break;
          pos = comma + 1;
        }
        if (tok.size() < 3) continue;
        // Skip header row: any non-numeric first field.
        if (tok[0].empty()) continue;
        bool numeric = true;
        for (size_t ci = 0; ci < tok[0].size(); ci++) {
          char c = tok[0][ci];
          if (!(c >= '0' && c <= '9') && c != '-' && c != '+') {
            numeric = false; break;
          }
        }
        if (!numeric) continue;
        BranchKey k;
        k.from = atoi(tok[0].c_str());
        k.to   = atoi(tok[1].c_str());
        k.ckt  = tok[2];
        monitorSet.insert(k);
      }
      if (world.rank() == 0) {
        printf("Monitor allowlist: %zu branches loaded from %s\n",
               monitorSet.size(), monitorBranchesFile.c_str());
      }
    }
  }
  // Endpoints of allowlisted branches; gates voltage rows under an allowlist.
  std::set<int> monitorBusSet;
  for (std::set<BranchKey>::const_iterator it = monitorSet.begin();
       it != monitorSet.end(); ++it) {
    monitorBusSet.insert(it->from);
    monitorBusSet.insert(it->to);
  }
  // Area/kV gate. Either-endpoint match for areas (catches tie-lines).
  // kV is gated on max(kv_from, kv_to) so a 138/13.8 stepdown counts as 138.
  // Empty area set / zero kV bound = unrestricted on that dimension.
  auto passesAreaKv = [&](int area_from, int area_to,
                          double kv_from, double kv_to) {
    if (!monitorAreas.empty()) {
      if (monitorAreas.find(area_from) == monitorAreas.end() &&
          monitorAreas.find(area_to)   == monitorAreas.end()) {
        return false;
      }
    }
    double kv_max = (kv_from > kv_to) ? kv_from : kv_to;
    if (monitorKvMin > 0.0 && kv_max < monitorKvMin) return false;
    if (monitorKvMax > 0.0 && kv_max > monitorKvMax) return false;
    return true;
  };
  // When monitorBranchesFile presents, it overrides area/kV criteria.
  bool haveAreaKvFilter = !monitorAreas.empty() ||
                          monitorKvMin > 0.0 ||
                          monitorKvMax > 0.0;
  if (!monitorSet.empty() && haveAreaKvFilter) {
    if (world.rank() == 0) {
      printf("WARNING: monitorBranchesFile is set; ignoring "
             "monitorAreas/monitorKvMin/monitorKvMax\n");
    }
    monitorAreas.clear();
    monitorKvMin = 0.0;
    monitorKvMax = 0.0;
    haveAreaKvFilter = false;
  }
  if (world.rank() == 0) {
    if (!monitorAreas.empty()) {
      printf("Monitor areas filter: %zu areas\n", monitorAreas.size());
    }
    if (monitorKvMin > 0.0 || monitorKvMax > 0.0) {
      printf("Monitor kV filter: min=%.2f max=%.2f (0 means unbounded)\n",
             monitorKvMin, monitorKvMax);
    }
  }
  auto busAreaKv = [&](int bus, int &area, double &kv) {
    std::map<int, BusAreaKv>::const_iterator it = bus_ak.find(bus);
    area = (it != bus_ak.end()) ? it->second.area   : 0;
    kv   = (it != bus_ak.end()) ? it->second.basekv : 0.0;
  };
  // Branch monitor predicate shared by every output path (ckt padding trimmed).
  auto branchMonitored = [&](int from, int to, const std::string &ckt) -> bool {
    if (!monitorSet.empty()) {
      BranchKey k;
      k.from = from; k.to = to; k.ckt = ckt;
      while (!k.ckt.empty() && (k.ckt[k.ckt.size()-1] == ' ' ||
                                k.ckt[k.ckt.size()-1] == '\t'))
        k.ckt.resize(k.ckt.size()-1);
      return monitorSet.find(k) != monitorSet.end();
    }
    if (!haveAreaKvFilter) return true;
    int af = 0, at = 0;
    double kf = 0.0, kt = 0.0;
    busAreaKv(from, af, kf);
    busAreaKv(to,   at, kt);
    return passesAreaKv(af, at, kf, kt);
  };
  // Bus counterpart: allowlist endpoint, else own area / base kV.
  auto busMonitored = [&](int bus) -> bool {
    if (!monitorSet.empty()) {
      return monitorBusSet.find(bus) != monitorBusSet.end();
    }
    if (!haveAreaKvFilter) return true;
    int area = 0;
    double kv = 0.0;
    busAreaKv(bus, area, kv);
    if (!monitorAreas.empty() && monitorAreas.find(area) == monitorAreas.end())
      return false;
    if (monitorKvMin > 0.0 && kv < monitorKvMin) return false;
    if (monitorKvMax > 0.0 && kv > monitorKvMax) return false;
    return true;
  };

  struct BaseFlow {
    double p_mw, q_mvar, mva, loading_pct;
    double base_rate, cont_rate;
    double v_from_pu, v_to_pu, ang_from_deg, ang_to_deg;
    double base_kv_from, base_kv_to;
    int    area_from, area_to;
  };
  std::map<BranchKey, BaseFlow> base_cache;
  std::string deltaPartPath;
  if (outputFormat == "csv_delta") {
    std::ostringstream oss;
    oss << outputFile << "_delta." << world.rank() << ".part";
    deltaPartPath = oss.str();
  }
  std::ofstream deltaPart;
  size_t deltaRowCount = 0;
  size_t deltaSkipCount = 0;

  // Convergence sidecar rows.
  struct ConvRow {
    int    event_idx;
    std::string name;
    std::string type;
    gridpack::utility::ConvergenceSummary cs;
    std::string status;
  };
  std::vector<ConvRow> localConvRows;
  // _convergence.csv: written for every outputFormat.
  bool emitConv = true;

  // Lambda: parse current solved flow_str/vr_str and stream one CSV row
  // per branch into the rank's .part file. Called once per converged case
  // (base + each contingency) on every task communicator; non-rank-0
  // task_comm members short-circuit after the collective.
  auto captureFlatRows = [&](int event_idx, const std::string &name,
                             bool emit, bool is_base) {
    std::vector<std::string> v_strs = pf_app.writeBusString("vr_str");
    std::vector<std::string> b_strs = pf_app.writeBranchString("flow_str");
    if (!emit || task_comm.rank() != 0) return;
    if (!flatPart.is_open()) {
      flatPart.open(flatPartPath.c_str(), std::ios::out | std::ios::trunc);
      flatPart << std::fixed;
    }
    std::map<int, std::pair<double,double> > vbymag_ang;
    for (size_t vi = 0; vi < v_strs.size(); vi++) {
      int    bus_id = 0, use_vmag = 0, changed = 0;
      double angle = 0.0, vmag = 0.0;
      if (sscanf(v_strs[vi].c_str(), "%d %lf %lf %d %d",
                 &bus_id, &angle, &vmag, &use_vmag, &changed) == 5) {
        vbymag_ang[bus_id] = std::make_pair(vmag, angle);
      }
    }
    char ct_name[24];
    std::strncpy(ct_name, name.c_str(), sizeof(ct_name) - 1);
    ct_name[sizeof(ct_name) - 1] = '\0';
    // Voltage PI and violations on monitored buses, contingency rows only.
    if (!is_base) {
      for (std::map<int,std::pair<double,double> >::const_iterator vit =
             vbymag_ang.begin(); vit != vbymag_ang.end(); ++vit) {
        if (!busMonitored(vit->first)) continue;
        double v_pu = vit->second.first;
        accumVoltagePi(ct_name, v_pu);
        if (v_pu <= 0.0 || !std::isfinite(v_pu)) continue;
        if (v_pu < Vmin || v_pu > Vmax) {
          emitVoltageViolation(event_idx, ct_name, vit->first, v_pu, Vmin, Vmax);
        }
      }
    }
    for (size_t bi = 0; bi < b_strs.size(); bi++) {
      char ckt_buf[16] = {0};
      int viol = 0;
      double p = 0.0, q = 0.0, perf = 0.0, ratea = 0.0;
      int from = 0, to = 0;
      if (sscanf(b_strs[bi].c_str(),
                 "%d %d %15s %lf %lf %lf %lf %d",
                 &from, &to, ckt_buf, &p, &q, &perf, &ratea, &viol) != 8) {
        continue;
      }
      char ckt[4];
      std::strncpy(ckt, ckt_buf, 3); ckt[3] = '\0';
      BranchKey mk;
      mk.from = from; mk.to = to; mk.ckt = ckt;
      while (!mk.ckt.empty() && mk.ckt[mk.ckt.size()-1] == ' ')
        mk.ckt.resize(mk.ckt.size()-1);
      if (!branchMonitored(from, to, mk.ckt)) continue;
      std::map<BranchKey, BranchRates>::const_iterator rIt = branch_rates.find(mk);
      double rate_sel = ratea;
      if (rIt != branch_rates.end()) {
        rate_sel = is_base ? rIt->second.rate_a : pickContRate(rIt->second);
      }
      double flow_mva    = std::sqrt(p*p + q*q);
      double loading_pct = (rate_sel > 0.0) ? (flow_mva / rate_sel) * 100.0 : 0.0;
      if (!is_base) accumBranchPi(ct_name, flow_mva, rate_sel);
      // Stream to _violations.csv for csv_flat runs (contingency rows only).
      if (!is_base && loading_pct > violationSeverityThreshold * 100.0) {
        double base_mva = 0.0;
        // No base_cache in csv_flat mode; look up in the persistent map built
        // by populateBaseCache-style capture below? We don't have one for
        // csv_flat, so base_mva stays 0 -- delta will just equal mva.
        emitBranchViolation(event_idx, ct_name, from, to, ckt,
                            flow_mva, rate_sel, loading_pct, base_mva);
      }
      std::map<int, std::pair<double,double> >::const_iterator vf =
        vbymag_ang.find(from);
      std::map<int, std::pair<double,double> >::const_iterator vt =
        vbymag_ang.find(to);
      double v_from       = (vf != vbymag_ang.end()) ? vf->second.first  : 0.0;
      double ang_from_deg = (vf != vbymag_ang.end()) ? vf->second.second : 0.0;
      double v_to         = (vt != vbymag_ang.end()) ? vt->second.first  : 0.0;
      double ang_to_deg   = (vt != vbymag_ang.end()) ? vt->second.second : 0.0;
      flatPart << event_idx << "," << ct_name << ","
               << from << "," << to << "," << ckt << ","
               << std::setprecision(4) << p << ","
               << std::setprecision(4) << q << ","
               << std::setprecision(4) << flow_mva << ","
               << std::setprecision(4) << rate_sel << ","
               << std::setprecision(2) << loading_pct << ","
               << viol << ","
               << std::setprecision(6) << v_from << ","
               << std::setprecision(6) << v_to << ","
               << std::setprecision(4) << ang_from_deg << ","
               << std::setprecision(4) << ang_to_deg
               << "\n";
      flatRowCount++;
    }
  };

  // Populate base_cache from current solved state. Called once after base
  // solve on every rank (csv_delta only); world.rank() == 0 is not special
  // here -- each rank caches the branches it sees on its task_comm so it
  // can join later in captureDeltaRows.
  auto populateBaseCache = [&]() {
    std::vector<std::string> v_strs = pf_app.writeBusString("vr_str");
    std::vector<std::string> b_strs = pf_app.writeBranchString("flow_str");
    if (task_comm.rank() != 0) return;
    std::map<int, std::pair<double,double> > vbymag_ang;
    for (size_t vi = 0; vi < v_strs.size(); vi++) {
      int    bus_id = 0, use_vmag = 0, changed = 0;
      double angle = 0.0, vmag = 0.0;
      if (sscanf(v_strs[vi].c_str(), "%d %lf %lf %d %d",
                 &bus_id, &angle, &vmag, &use_vmag, &changed) == 5) {
        vbymag_ang[bus_id] = std::make_pair(vmag, angle);
      }
    }
    for (size_t bi = 0; bi < b_strs.size(); bi++) {
      char ckt_buf[16] = {0};
      int viol = 0;
      double p = 0.0, q = 0.0, perf = 0.0, ratea = 0.0;
      int from = 0, to = 0;
      if (sscanf(b_strs[bi].c_str(),
                 "%d %d %15s %lf %lf %lf %lf %d",
                 &from, &to, ckt_buf, &p, &q, &perf, &ratea, &viol) != 8) {
        continue;
      }
      BranchKey k;
      k.from = from; k.to = to;
      k.ckt  = std::string(ckt_buf);
      // Strip trailing spaces from ckt so the key matches what flow_str
      // returns later (sscanf %15s already trims leading whitespace).
      while (!k.ckt.empty() && k.ckt[k.ckt.size()-1] == ' ') k.ckt.resize(k.ckt.size()-1);
      if (!branchMonitored(from, to, k.ckt)) continue;
      double base_rate = ratea, cont_rate = ratea;
      std::map<BranchKey, BranchRates>::const_iterator rIt = branch_rates.find(k);
      if (rIt != branch_rates.end()) {
        base_rate = rIt->second.rate_a;
        cont_rate = pickContRate(rIt->second);
      }
      BaseFlow bf;
      bf.p_mw        = p;
      bf.q_mvar      = q;
      bf.mva         = std::sqrt(p*p + q*q);
      bf.base_rate   = base_rate;
      bf.cont_rate   = cont_rate;
      bf.loading_pct = (base_rate > 0.0) ? (bf.mva / base_rate) * 100.0 : 0.0;
      std::map<int, std::pair<double,double> >::const_iterator vf =
        vbymag_ang.find(from);
      std::map<int, std::pair<double,double> >::const_iterator vt =
        vbymag_ang.find(to);
      bf.v_from_pu     = (vf != vbymag_ang.end()) ? vf->second.first  : 0.0;
      bf.ang_from_deg  = (vf != vbymag_ang.end()) ? vf->second.second : 0.0;
      bf.v_to_pu       = (vt != vbymag_ang.end()) ? vt->second.first  : 0.0;
      bf.ang_to_deg    = (vt != vbymag_ang.end()) ? vt->second.second : 0.0;
      busAreaKv(from, bf.area_from, bf.base_kv_from);
      busAreaKv(to,   bf.area_to,   bf.base_kv_to);
      base_cache[k] = bf;
    }
  };

  // Wide-form (base+cont on same row) capture for csv_delta. Mirrors
  // captureFlatRows but joins each branch with base_cache. Branches not
  // in base_cache are counted in deltaSkipCount and skipped silently.
  auto captureDeltaRows = [&](int event_idx,
                              const gridpack::powerflow::Contingency &evt,
                              bool emit) {
    std::vector<std::string> v_strs = pf_app.writeBusString("vr_str");
    std::vector<std::string> b_strs = pf_app.writeBranchString("flow_str");
    if (!emit || task_comm.rank() != 0) return;
    if (!deltaPart.is_open()) {
      deltaPart.open(deltaPartPath.c_str(), std::ios::out | std::ios::trunc);
      deltaPart << std::fixed;
    }
    // cont_event_facility: disabled -- join event_idx against
    // <outputFile>_contingencies.csv instead, which names every outaged element
    // rather than the first plus "(+N more)". Kept commented in case the column
    // is wanted back; uncomment this block, the emission below and the header
    // field together. The [area] lookup needs a complete bus_meta (groupSize=1).
    // std::string facility;
    // // clean2Char pads ids to two chars; trim so the label has no stray space.
    // auto rtrimId = [](const std::string &in) -> std::string {
    //   std::string t = in;
    //   while (!t.empty() && (t[t.size()-1] == ' ' || t[t.size()-1] == '\t'))
    //     t.resize(t.size()-1);
    //   return t;
    // };
    // // "[area] <element ids>" for both kinds; the type column says which.
    // if (evt.p_type == Branch && !evt.p_from.empty()) {
    //   int outFrom = evt.p_from[0];
    //   int area = 0;
    //   std::map<int, BusMeta>::const_iterator mf = bus_meta.find(outFrom);
    //   if (mf != bus_meta.end()) area = mf->second.area;
    //   char buf[64];
    //   snprintf(buf, sizeof(buf), "[%d] %d %d %s",
    //            area, outFrom, evt.p_to[0], rtrimId(evt.p_ckt[0]).c_str());
    //   facility = buf;
    //   if (evt.p_from.size() > 1) {
    //     char suf[24];
    //     snprintf(suf, sizeof(suf), " (+%zu more)", evt.p_from.size() - 1);
    //     facility += suf;
    //   }
    // } else if (evt.p_type == Generator && !evt.p_busid.empty()) {
    //   int outBus = evt.p_busid[0];
    //   int area = 0;
    //   std::map<int, BusMeta>::const_iterator mg = bus_meta.find(outBus);
    //   if (mg != bus_meta.end()) area = mg->second.area;
    //   char buf[64];
    //   snprintf(buf, sizeof(buf), "[%d] %d %s",
    //            area, outBus, rtrimId(evt.p_genid[0]).c_str());
    //   facility = buf;
    //   if (evt.p_busid.size() > 1) {
    //     char suf[24];
    //     snprintf(suf, sizeof(suf), " (+%zu more)", evt.p_busid.size() - 1);
    //     facility += suf;
    //   }
    // }
    std::string ct_name = evt.p_name;
    while (!ct_name.empty() && ct_name[ct_name.size()-1] == ' ')
      ct_name.resize(ct_name.size()-1);
    const char *type_str = (evt.p_type == Branch) ? "branch" : "generator";
    std::map<int, std::pair<double,double> > vbymag_ang;
    for (size_t vi = 0; vi < v_strs.size(); vi++) {
      int    bus_id = 0, use_vmag = 0, changed = 0;
      double angle = 0.0, vmag = 0.0;
      if (sscanf(v_strs[vi].c_str(), "%d %lf %lf %d %d",
                 &bus_id, &angle, &vmag, &use_vmag, &changed) == 5) {
        vbymag_ang[bus_id] = std::make_pair(vmag, angle);
      }
    }
    // Voltage PI and violations on monitored buses.
    for (std::map<int,std::pair<double,double> >::const_iterator vit =
           vbymag_ang.begin(); vit != vbymag_ang.end(); ++vit) {
      if (!busMonitored(vit->first)) continue;
      double v_pu = vit->second.first;
      accumVoltagePi(ct_name, v_pu);
      if (v_pu <= 0.0 || !std::isfinite(v_pu)) continue;
      if (v_pu < Vmin || v_pu > Vmax) {
        emitVoltageViolation(event_idx, ct_name, vit->first, v_pu, Vmin, Vmax);
      }
    }
    for (size_t bi = 0; bi < b_strs.size(); bi++) {
      char ckt_buf[16] = {0};
      int viol = 0;
      double p = 0.0, q = 0.0, perf = 0.0, ratea = 0.0;
      int from = 0, to = 0;
      if (sscanf(b_strs[bi].c_str(),
                 "%d %d %15s %lf %lf %lf %lf %d",
                 &from, &to, ckt_buf, &p, &q, &perf, &ratea, &viol) != 8) {
        continue;
      }
      BranchKey k;
      k.from = from; k.to = to;
      k.ckt  = std::string(ckt_buf);
      while (!k.ckt.empty() && k.ckt[k.ckt.size()-1] == ' ')
        k.ckt.resize(k.ckt.size()-1);
      if (!branchMonitored(from, to, k.ckt)) continue;
      std::map<BranchKey, BaseFlow>::const_iterator it = base_cache.find(k);
      if (it == base_cache.end()) { deltaSkipCount++; continue; }
      const BaseFlow &bf = it->second;
      double cont_mva     = std::sqrt(p*p + q*q);
      double cont_loading = (bf.cont_rate > 0.0) ? (cont_mva / bf.cont_rate) * 100.0 : 0.0;
      accumBranchPi(ct_name, cont_mva, bf.cont_rate);
      // Stream to _violations.csv (delta path knows base_mva already).
      if (cont_loading > violationSeverityThreshold * 100.0) {
        emitBranchViolation(event_idx, ct_name, from, to, k.ckt,
                            cont_mva, bf.cont_rate, cont_loading, bf.mva);
      }
      std::map<int, std::pair<double,double> >::const_iterator vf =
        vbymag_ang.find(from);
      std::map<int, std::pair<double,double> >::const_iterator vt =
        vbymag_ang.find(to);
      double v_from_c = (vf != vbymag_ang.end()) ? vf->second.first  : 0.0;
      double a_from_c = (vf != vbymag_ang.end()) ? vf->second.second : 0.0;
      double v_to_c   = (vt != vbymag_ang.end()) ? vt->second.first  : 0.0;
      double a_to_c   = (vt != vbymag_ang.end()) ? vt->second.second : 0.0;
      double d_ang_b  = bf.ang_from_deg - bf.ang_to_deg;
      double d_ang_c  = a_from_c - a_to_c;
      // Across-branch drop, same convention as the angle deltas above.
      double d_v_b    = bf.v_from_pu - bf.v_to_pu;
      double d_v_c    = v_from_c - v_to_c;
      deltaPart << event_idx << "," << ct_name << "," << type_str << ","
                << from << "," << to << "," << k.ckt << ","
                << std::setprecision(2) << bf.base_kv_from << ","
                << std::setprecision(2) << bf.base_kv_to   << ","
                << bf.area_from << "," << bf.area_to << ","
                << std::setprecision(4) << bf.base_rate << ","
                << std::setprecision(4) << bf.cont_rate << ","
                << std::setprecision(4) << bf.p_mw   << ","
                << std::setprecision(4) << p         << ","
                << std::setprecision(4) << bf.q_mvar << ","
                << std::setprecision(4) << q         << ","
                << std::setprecision(4) << bf.mva    << ","
                << std::setprecision(4) << cont_mva  << ","
                << std::setprecision(2) << bf.loading_pct << ","
                << std::setprecision(2) << cont_loading  << ","
                << std::setprecision(6) << bf.v_from_pu << ","
                << std::setprecision(6) << v_from_c     << ","
                << std::setprecision(6) << bf.v_to_pu   << ","
                << std::setprecision(6) << v_to_c       << ","
                << std::setprecision(4) << bf.ang_from_deg << ","
                << std::setprecision(4) << a_from_c        << ","
                << std::setprecision(4) << bf.ang_to_deg   << ","
                << std::setprecision(4) << a_to_c          << ","
                << std::setprecision(6) << d_v_b   << ","
                << std::setprecision(6) << d_v_c   << ","
                << std::setprecision(4) << d_ang_b << ","
                << std::setprecision(4) << d_ang_c
             // << "," << facility            // cont_event_facility: disabled
                << "\n";
      deltaRowCount++;
    }
  };

  //  Set minimum and maximum voltage limits on all buses
  pf_app.setVoltageLimits(Vmin, Vmax);
  // Route CA violation checks and loadingPercent through the same rating tier.
  pf_app.setContingencyRating(contingencyRating);
  // Solve the base power flow on every task communicator. Abort if it fails.
  bool baseSolveOk = false;
  try {
    baseSolveOk = pf_app.solve();
    if (baseSolveOk && check_Qlim && !pf_app.checkQlimViolations()) {
      baseSolveOk = pf_app.solve();
    }
  } catch (const std::exception &e) {
    if (world.rank() == 0) {
      printf("ERROR: base-case solve threw exception: %s\n", e.what());
    }
    baseSolveOk = false;
  } catch (...) {
    if (world.rank() == 0) {
      printf("ERROR: base-case solve threw unknown exception\n");
    }
    baseSolveOk = false;
  }
  if (!baseSolveOk) {
    if (world.rank() == 0) {
      gridpack::utility::ConvergenceSummary cs = pf_app.getConvergence();
      printf("ERROR: base case did not converge "
             "(iterations=%d, final_tol=%.6e, "
             "max_p_bus=%d max_p_mismatch=%.4f, "
             "max_q_bus=%d max_q_mismatch=%.4f). "
             "Aborting contingency analysis.\n",
             cs.iterations, cs.finalTolerance,
             cs.finalMismatch.maxPBus, cs.finalMismatch.maxPMismatch,
             cs.finalMismatch.maxQBus, cs.finalMismatch.maxQMismatch);
    }
    world.barrier();
    MPI_Abort(static_cast<MPI_Comm>(world), 1);
  }
  // Suppress voltage violations already present at base.
  pf_app.ignoreVoltageViolations();

  // Flag non-monitored elements "ignore" so the violation checks and .out
  // listings follow the filter.
  if (haveMonitorFilter) {
    long cnt[4] = { 0, 0, 0, 0 };   // mon buses, buses, mon elems, elems
    int nBus = pf_network->numBuses();
    for (int i = 0; i < nBus; i++) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>(pf_network->getBus(i).get());
      if (!bus) continue;
      bool mon = busMonitored(pf_network->getOriginalBusIndex(i));
      if (!mon) bus->setIgnore(true);
      if (pf_network->getActiveBus(i)) { cnt[1]++; if (mon) cnt[0]++; }
    }
    int nBranch = pf_network->numBranches();
    for (int i = 0; i < nBranch; i++) {
      gridpack::powerflow::PFBranch *br =
        dynamic_cast<gridpack::powerflow::PFBranch*>(pf_network->getBranch(i).get());
      if (!br) continue;
      int from = br->getBus1OriginalIndex();
      int to   = br->getBus2OriginalIndex();
      std::vector<std::string> tags = br->getLineTags();
      bool active = pf_network->getActiveBranch(i);
      for (size_t t = 0; t < tags.size(); t++) {
        bool mon = branchMonitored(from, to, tags[t]);
        if (!mon) br->setIgnore(tags[t], true);
        if (active) { cnt[3]++; if (mon) cnt[2]++; }
      }
    }
    long tot[4] = { 0, 0, 0, 0 };
    MPI_Allreduce(cnt, tot, 4, MPI_LONG, MPI_SUM,
                  static_cast<MPI_Comm>(task_comm));
    if (world.rank() == 0) {
      printf("Monitor filter: %ld of %ld branch elements and %ld of %ld buses "
             "monitored\n", tot[2], tot[3], tot[0], tot[1]);
      if (tot[2] == 0 && tot[0] == 0) {
        printf("WARNING: monitor filter matches nothing; all outputs will be "
               "empty. Check monitorAreas/monitorKvMin/monitorKvMax/"
               "monitorBranchesFile against the case.\n");
      }
    }
  }
  // Drop non-monitored elements from a collected result set.
  auto filterResults = [&](gridpack::utility::PowerFlowResults &r) {
    if (!haveMonitorFilter) return;
    std::vector<gridpack::utility::BusResult> buses;
    for (size_t i = 0; i < r.buses.size(); i++) {
      if (busMonitored(r.buses[i].busId)) buses.push_back(r.buses[i]);
    }
    r.buses.swap(buses);
    std::vector<gridpack::utility::BranchResult> branches;
    for (size_t i = 0; i < r.branches.size(); i++) {
      const gridpack::utility::BranchResult &b = r.branches[i];
      if (branchMonitored(b.fromBus, b.toBus, b.circuitId)) branches.push_back(b);
    }
    r.branches.swap(branches);
    std::vector<gridpack::utility::GeneratorResult> gens;
    for (size_t i = 0; i < r.generators.size(); i++) {
      if (busMonitored(r.generators[i].busId)) gens.push_back(r.generators[i]);
    }
    r.generators.swap(gens);
  };

  // Collect base case results for export. csv_flat captures rows directly
  // in the hot loop and skips the heavyweight collectResults() path.
  gridpack::utility::PowerFlowResults baseCaseResults;
  if (outputFormat == "json" || outputFormat == "csv" ||
      outputFormat == "text") {
    baseCaseResults = pf_app.collectResults();
    filterResults(baseCaseResults);
  }
  if (outputFormat == "csv_flat") {
    // The base case is replicated on every task communicator. captureFlatRows
    // calls writeBusString/writeBranchString which are task_comm collectives,
    // so every task_comm participates -- but only world rank 0 emits rows so
    // the base case isn't duplicated in the final file.
    captureFlatRows(0, std::string("base_case"), world.rank() == 0, true);
  }
  if (outputFormat == "csv_delta") {
    // Cache base-case branch state on every rank for the contingency join.
    populateBaseCache();
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

  // event_idx lookup table: one row per contingency (0 = base case); N-k
  // element ids are ';'-separated within the columns.
  if (world.rank() == 0) {
    std::string ctgFile = outputFile + "_contingencies.csv";
    std::ofstream cout_ctg(ctgFile.c_str(), std::ios::out | std::ios::trunc);
    cout_ctg << "event_idx,contingency,type,n_elements,"
                "from_bus,to_bus,circuit_id,gen_bus,gen_id\n";
    size_t ctgRows = 0;
    cout_ctg << "0,base_case,base,0,,,,,\n";
    ctgRows++;
    // Trim clean2Char padding so ids join against the other CSVs.
    auto rtrim = [](const std::string &in) -> std::string {
      std::string t = in;
      while (!t.empty() && (t[t.size()-1] == ' ' || t[t.size()-1] == '\t'))
        t.resize(t.size()-1);
      return t;
    };
    for (size_t ei = 0; ei < events.size(); ei++) {
      const gridpack::powerflow::Contingency &e = events[ei];
      int event_idx = static_cast<int>(ei) + 1;
      std::string nm = rtrim(e.p_name);
      size_t n = 0;
      const char *ty = "unknown";
      if (e.p_type == Branch) { n = e.p_from.size(); ty = "branch"; }
      else if (e.p_type == Generator) { n = e.p_busid.size(); ty = "generator"; }
      std::ostringstream c_from, c_to, c_ckt, c_gbus, c_gid;
      for (size_t j = 0; j < n; j++) {
        const char *sep = (j > 0) ? ";" : "";
        if (e.p_type == Branch) {
          c_from << sep << e.p_from[j];
          c_to   << sep << e.p_to[j];
          c_ckt  << sep << rtrim(e.p_ckt[j]);
        } else {
          c_gbus << sep << e.p_busid[j];
          c_gid  << sep << rtrim(e.p_genid[j]);
        }
      }
      // Empty events still get a row so every event_idx decodes.
      cout_ctg << event_idx << "," << nm << "," << ty << "," << n << ","
               << c_from.str() << "," << c_to.str() << "," << c_ckt.str() << ","
               << c_gbus.str() << "," << c_gid.str() << "\n";
      ctgRows++;
    }
    cout_ctg.close();
    printf("[contingencies] wrote %zu rows to %s\n", ctgRows, ctgFile.c_str());
  }

  // Print contingency details (gated on printCalcFiles; noisy for large lists)
  if (print_calcs && world.rank() == 0) {
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
  // StatBlock objects and the per-case scratch vectors live across the
  // contingency loop so they are declared up here, regardless of whether
  // statistics output is enabled.
  boost::scoped_ptr<gridpack::analysis::StatBlock> vmag_stats;
  boost::scoped_ptr<gridpack::analysis::StatBlock> vang_stats;
  boost::scoped_ptr<gridpack::analysis::StatBlock> pgen_stats;
  boost::scoped_ptr<gridpack::analysis::StatBlock> qgen_stats;
  boost::scoped_ptr<gridpack::analysis::StatBlock> pflow_stats;
  boost::scoped_ptr<gridpack::analysis::StatBlock> qflow_stats;
  boost::scoped_ptr<gridpack::analysis::StatBlock> perf_stats;
  std::vector<std::string> v_vals;
  int nsize = 0;
  std::vector<double> vmag, vang, pgen, qgen, pflow, qflow, perf;
  std::vector<int> mask, mag_mask;
  int t_store = timer->createCategory("Store Statistics");
  if (write_stats) {
    timer->start(t_store);
    v_vals = pf_app.writeBusString("vr_str");
    nsize = v_vals.size();
    std::vector<int> mag_ids;
    std::vector<int> ids;
    std::vector<std::string> mag_tags;
    std::vector<std::string> tags;
    // Find bus IDs and create a dummy tag label and get voltage magnitude
    // and angle for base case
    for (i=0; i<nsize; i++) {
      std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
      if (!busMonitored(atoi(tokens[0].c_str()))) continue;
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
    // Create StatBlock objects for voltage magnitude and angles and add
    // bus IDs to it as well as base case values
    vmag_stats.reset(new gridpack::analysis::StatBlock(world,nmags,ntasks+1));
    vang_stats.reset(new gridpack::analysis::StatBlock(world,nbus,ntasks+1));
    if (world.rank() == 0) {
      vmag_stats->addRowLabels(mag_ids, mag_tags);
      vang_stats->addRowLabels(ids, tags);
      vmag_stats->addColumnValues(0,vmag,mag_mask);
      vang_stats->addColumnValues(0,vang,mask);
    }
    // Get generator power information
    v_vals.clear();
    ids.clear();
    tags.clear();
    mask.clear();
    v_vals = pf_app.writeBusString("power");
    nsize = v_vals.size();
    // Find bus IDs and tags for generators and evaluate Pg and Qg for base case
    for (i=0; i<nsize; i++) {
      std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
      if (tokens.size()%4 != 0) {
        printf("Incorrect generator listing\n");
        continue;
      }
      int ngen = tokens.size()/4;
      for (j=0; j<ngen; j++) {
        if (!busMonitored(atoi(tokens[j*4].c_str()))) continue;
        ids.push_back(atoi(tokens[j*4].c_str()));
        tags.push_back(tokens[j*4+1]);
        pgen.push_back(atof(tokens[j*4+2].c_str()));
        qgen.push_back(atof(tokens[j*4+3].c_str()));
        mask.push_back(1);
      }
    }
    nsize = pgen.size();
    world.max(&nsize,1);
    // Create StatBlock objects for Pg and Qg and add labels and base case values
    pgen_stats.reset(new gridpack::analysis::StatBlock(world,nsize,ntasks+1));
    qgen_stats.reset(new gridpack::analysis::StatBlock(world,nsize,ntasks+1));
    if (world.rank() == 0) {
      pgen_stats->addRowLabels(ids, tags);
      qgen_stats->addRowLabels(ids, tags);
      pgen_stats->addColumnValues(0,pgen,mask);
      qgen_stats->addColumnValues(0,qgen,mask);
    }

    // Find flow parameters for all branch lines
    v_vals.clear();
    ids.clear();
    tags.clear();
    mask.clear();
    std::vector<int> id1;
    std::vector<int> id2;
    std::vector<double> pmin, pmax;
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
        if (!branchMonitored(atoi(tokens[j*8].c_str()),
                             atoi(tokens[j*8+1].c_str()),
                             tokens[j*8+2])) continue;
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
    // Create StatBlock objects for flow parameters and add labels and base case
    // values
    pflow_stats.reset(new gridpack::analysis::StatBlock(world,nsize,ntasks+1));
    qflow_stats.reset(new gridpack::analysis::StatBlock(world,nsize,ntasks+1));
    perf_stats.reset(new gridpack::analysis::StatBlock(world,nsize,ntasks+1));
    if (world.rank() == 0) {
      pflow_stats->addRowLabels(id1, id2, tags);
      qflow_stats->addRowLabels(id1, id2, tags);
      perf_stats->addRowLabels(id1, id2, tags);
      pflow_stats->addColumnValues(0,pflow,mask);
      qflow_stats->addColumnValues(0,qflow,mask);
      perf_stats->addColumnValues(0,perf,mask);
      pflow_stats->addRowMinValue(pmin);
      qflow_stats->addRowMinValue(pmin);
      pflow_stats->addRowMaxValue(pmax);
      qflow_stats->addRowMaxValue(pmax);
    }
    timer->stop(t_store);
  }
  if (check_Qlim) pf_app.clearQlimViolations();
  // Clear any Q limit warnings from base case before starting contingencies
  gridpack::powerflow::PFBus::clearQlimWarnings();

  // Local contingency results storage for JSON/CSV export
  std::vector<gridpack::utility::ContingencyResult> localContingencies;

  // Base-case per-element MVA lookup, keyed on (from,to,ckt).
  // Populated on rank 0 of each task_comm (only place collectResults ran).
  // Used to fill BranchViolation.baseMva/deltaMva during contingency reporting.
  std::map<BranchKey, double> baseMvaByKey;
  if (outputFormat == "json" || outputFormat == "csv" ||
      outputFormat == "text") {
    for (size_t bi = 0; bi < baseCaseResults.branches.size(); bi++) {
      const gridpack::utility::BranchResult &br = baseCaseResults.branches[bi];
      BranchKey k; k.from = br.fromBus; k.to = br.toBus; k.ckt = br.circuitId;
      while (!k.ckt.empty() && k.ckt[k.ckt.size()-1] == ' ') k.ckt.resize(k.ckt.size()-1);
      double mva = (br.mvaFrom > br.mvaTo) ? br.mvaFrom : br.mvaTo;
      baseMvaByKey[k] = mva;
    }
  }

  // Fill BranchViolation/VoltageViolation arrays from a solved ct result,
  // and stream the same rows to <outputFile>_violations.<rank>.part.
  // event_idx of 0 is reserved for base case; contingencies get task_id+1.
  auto populateViolations = [&](gridpack::utility::ContingencyResult &ct,
                                int event_idx) {
    const double threshPct = violationSeverityThreshold * 100.0;
    // Gate PI on task_comm rank 0 so groupSize>1 doesn't double-count.
    bool accumPi = (task_comm.rank() == 0);
    for (size_t bi = 0; bi < ct.solution.branches.size(); bi++) {
      const gridpack::utility::BranchResult &br = ct.solution.branches[bi];
      if (accumPi) {
        double mva_br = (br.mvaFrom > br.mvaTo) ? br.mvaFrom : br.mvaTo;
        accumBranchPi(ct.name, mva_br, br.rateSelected);
      }
      if (br.loadingPercent <= threshPct) continue;
      gridpack::utility::BranchViolation v;
      v.fromBus = br.fromBus;
      v.toBus = br.toBus;
      v.circuitId = br.circuitId;
      v.mva = (br.mvaFrom > br.mvaTo) ? br.mvaFrom : br.mvaTo;
      v.rate = br.rateSelected;
      v.loadingPercent = br.loadingPercent;
      BranchKey k; k.from = br.fromBus; k.to = br.toBus; k.ckt = br.circuitId;
      while (!k.ckt.empty() && k.ckt[k.ckt.size()-1] == ' ') k.ckt.resize(k.ckt.size()-1);
      std::map<BranchKey, double>::const_iterator it = baseMvaByKey.find(k);
      v.baseMva = (it != baseMvaByKey.end()) ? it->second : 0.0;
      v.deltaMva = v.mva - v.baseMva;
      v.severity = (br.loadingPercent >= 105.0) ? "critical" : "warning";
      ct.branchViolations.push_back(v);
      emitBranchViolation(event_idx, ct.name, v.fromBus, v.toBus, v.circuitId,
                          v.mva, v.rate, v.loadingPercent, v.baseMva);
    }
    for (size_t bi = 0; bi < ct.solution.buses.size(); bi++) {
      const gridpack::utility::BusResult &b = ct.solution.buses[bi];
      double v_pu = b.voltage;
      if (v_pu <= 0.0) continue;   // Skip isolated / not-solved buses.
      // Voltage PI accrues on every energized bus (textbook form),
      // gated on task_comm rank 0 to avoid double-count under groupSize>1.
      if (accumPi) accumVoltagePi(ct.name, v_pu);
      bool lo = v_pu < Vmin, hi = v_pu > Vmax;
      if (!lo && !hi) continue;
      gridpack::utility::VoltageViolation vv;
      vv.busId = b.busId;
      vv.vPu = v_pu;
      vv.limitLow = Vmin;
      vv.limitHigh = Vmax;
      vv.deviationPu = lo ? (v_pu - Vmin) : (v_pu - Vmax);
      double dev = std::abs(vv.deviationPu);
      vv.severity = (dev >= 0.05) ? "critical" : "warning";
      ct.voltageViolations.push_back(vv);
      emitVoltageViolation(event_idx, ct.name, vv.busId, vv.vPu,
                           vv.limitLow, vv.limitHigh);
    }
  };

  // Convergence row recorder; indexes events[task_id].
  auto recordConv = [&](int task_id, const char *status,
                        const std::string &) {
    if (!emitConv) return;
    if (task_comm.rank() != 0) return;
    ConvRow r;
    r.event_idx = task_id + 1;
    r.name      = events[task_id].p_name;
    r.type      = (events[task_id].p_type == Branch) ? "branch" : "generator";
    r.cs        = pf_app.getConvergence();
    r.status    = status;
    localConvRows.push_back(r);
  };

  // Evaluate contingencies using the task manager
  int task_id;
  char sbuf[512];
  // nextTask returns the same task_id on all processors in task_comm. When the
  // calculation runs out of task, nextTask will return false.
  while (taskmgr.nextTask(task_comm, &task_id)) {
    if (print_calcs) printf("Executing task %d on process %d\n",task_id,world.rank());
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
        if (print_calcs) printf("p[%d] Line: (from) %d (to) %d (line) \'%s\'\n",
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
        if (print_calcs) printf("p[%d] Generator: (bus) %d (generator ID) \'%s\'\n",
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
    // Skip power flow if contingency setup failed (no valid slack) or islanding detected
    bool slackCapacityOk = true;  // Will be checked after solve
    bool solveOk = false;
    if (contingencyFound && !islandDetected) {
      try {
        solveOk = pf_app.solve();
        if (solveOk && check_Qlim && !pf_app.checkQlimViolations()) {
          pf_app.solve();
        }
      } catch (const std::exception& e) {
        printf("p[%d] hit exception: %s\n", world.rank(), e.what());
        printf("Solver failure\n");
        solveOk = false;
      } catch (...) {
        printf("p[%d] hit unknown exception during solve\n", world.rank());
        solveOk = false;
      }
    }
    if (solveOk) {
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
        if (outputFormat == "json" || outputFormat == "csv") {
          gridpack::utility::ContingencyResult ctResult;
          ctResult.name = events[task_id].p_name;
          ctResult.type = (events[task_id].p_type == Branch) ? "branch" : "generator";
          ctResult.hasVoltageViolation = false;
          ctResult.hasBranchViolation = false;
          ctResult.solution.convergence = pf_app.getConvergence();
          ctResult.solution.convergence.converged = false;
          localContingencies.push_back(ctResult);
        }
        recordConv(task_id, "SLACK_OVERLOAD", std::string());
        sprintf(sbuf,"\nInsufficient generation capacity for contingency %s\n",
            events[task_id].p_name.c_str());
        if (print_calcs) pf_app.print(sbuf);
      } else {
        // Power flow solved and slack within capacity
        // If power flow solution is successful, write out voltages and currents
        if (print_calcs) pf_app.write();
        // Check for violations
        bool ok1 = pf_app.checkVoltageViolations();
        bool ok2 = pf_app.checkLineOverloadViolations();
        bool ok = ok1 && ok2;
        // text mode runs the summary path but discards the per-ct struct.
        if (outputFormat == "json" || outputFormat == "csv" ||
            outputFormat == "text") {
          gridpack::utility::ContingencyResult ctResult;
          ctResult.name = events[task_id].p_name;
          ctResult.type = (events[task_id].p_type == Branch) ? "branch" : "generator";
          ctResult.hasVoltageViolation = !ok1;
          ctResult.hasBranchViolation = !ok2;
          ctResult.solution = pf_app.collectResults();
          filterResults(ctResult.solution);
          populateViolations(ctResult, static_cast<int>(task_id) + 1);
          if (outputFormat != "text") localContingencies.push_back(ctResult);
        }
        if (outputFormat == "csv_flat") {
          captureFlatRows(task_id + 1, events[task_id].p_name, true, false);
        }
        if (outputFormat == "csv_delta") {
          captureDeltaRows(task_id + 1, events[task_id], true);
        }
        recordConv(task_id, "OK", std::string());
      // Include results of violation checks in output
      if (ok) {
        sprintf(sbuf,"\nNo violation for contingency %s\n",
            events[task_id].p_name.c_str());
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
        // Keep in step with the row format in PFBranch::serialWrite("flow").
        sprintf(sbuf,"\nBranch Violation for contingency %s\n"
            "  From Bus    To Bus   CKT       P_from       Q_from"
            "     MVA_from         P_to         Q_to       MVA_to"
            "       Rate   Loading%%\n",
            events[task_id].p_name.c_str());
      } else if (!ok) {
        sprintf(sbuf,"\nNo Branch Violation for contingency %s\n",
            events[task_id].p_name.c_str());
      }

      if (print_calcs) pf_app.print(sbuf);
      if (print_calcs) pf_app.writeCABranch();
      // Get strings of data from power flow calculation and parse them to
      // extract numerical values. Store these values in vectors and then
      // add them to StatBlock objects
      if (write_stats) {
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
          if (!busMonitored(atoi(tokens[0].c_str()))) continue;
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
        if (task_comm.rank() == 0) {
          vmag_stats->addColumnValues(task_id+1,vmag,mag_mask);
          vang_stats->addColumnValues(task_id+1,vang,mask);
        }
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
            if (!busMonitored(atoi(tokens[j*4].c_str()))) continue;
            pgen.push_back(atof(tokens[j*4+2].c_str()));
            qgen.push_back(atof(tokens[j*4+3].c_str()));
            mask.push_back(1);
          }
        }
        if (task_comm.rank() == 0) {
          pgen_stats->addColumnValues(task_id+1,pgen,mask);
          qgen_stats->addColumnValues(task_id+1,qgen,mask);
        }
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
            if (!branchMonitored(atoi(tokens[j*8].c_str()),
                                 atoi(tokens[j*8+1].c_str()),
                                 tokens[j*8+2])) continue;
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
        if (task_comm.rank() == 0) {
          pflow_stats->addColumnValues(task_id+1,pflow,mask);
          qflow_stats->addColumnValues(task_id+1,qflow,mask);
          perf_stats->addColumnValues(task_id+1,perf,mask);
        }
        timer->stop(t_store);
      }
        // Note: clearQlimViolations() moved after unSetContingency() below
      }  // end slackCapacityOk block
    } else {
      if (outputFormat == "json" || outputFormat == "csv") {
        gridpack::utility::ContingencyResult ctResult;
        ctResult.name = events[task_id].p_name;
        ctResult.type = (events[task_id].p_type == Branch) ? "branch" : "generator";
        ctResult.hasVoltageViolation = false;
        ctResult.hasBranchViolation = false;
        ctResult.solution.convergence = pf_app.getConvergence();
        ctResult.solution.convergence.converged = false;
        localContingencies.push_back(ctResult);
      }
      {
        const char *st;
        if (islandDetected) {
          st = "ISLANDED";
        } else if (!contingencyFound) {
          st = "NO_SLACK";
        } else {
          st = "DIVERGED";
        }
        recordConv(task_id, st, std::string());
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
      if (write_stats) {
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
          if (!busMonitored(atoi(tokens[0].c_str()))) continue;
          int not_isolated = atoi(tokens[3].c_str());
          if (not_isolated == 1) {
            vmag.push_back(0.0);
            mag_mask.push_back(0);
          }
          vang.push_back(0.0);
          mask.push_back(0);
        }
        if (task_comm.rank() == 0) {
          vmag_stats->addColumnValues(task_id+1,vmag,mag_mask);
          vang_stats->addColumnValues(task_id+1,vang,mask);
        }
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
            if (!busMonitored(atoi(tokens[j*4].c_str()))) continue;
            pgen.push_back(0.0);
            qgen.push_back(0.0);
            mask.push_back(0);
          }
        }
        if (task_comm.rank() == 0) {
          pgen_stats->addColumnValues(task_id+1,pgen,mask);
          qgen_stats->addColumnValues(task_id+1,qgen,mask);
        }
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
            if (!branchMonitored(atoi(tokens[j*8].c_str()),
                                 atoi(tokens[j*8+1].c_str()),
                                 tokens[j*8+2])) continue;
            pflow.push_back(0.0);
            qflow.push_back(0.0);
            perf.push_back(0.0);
            mask.push_back(0);
          }
        }
        if (task_comm.rank() == 0) {
          pflow_stats->addColumnValues(task_id+1,pflow,mask);
          qflow_stats->addColumnValues(task_id+1,qflow,mask);
          perf_stats->addColumnValues(task_id+1,perf,mask);
        }
        timer->stop(t_store);
      }
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
  // csv_flat / csv_delta: each rank streamed rows to its .part file during
  // the loop. Close, sync, then world rank 0 writes header + concatenates.
  if (outputFormat == "csv_flat") {
    if (flatPart.is_open()) flatPart.close();
  }
  if (outputFormat == "csv_delta") {
    if (deltaPart.is_open()) deltaPart.close();
  }
  if (wantBusSidecar) {
    world.sync();
    if (world.rank() == 0) {
      const size_t BUFSZ = 1 << 20;
      std::vector<char> buf(BUFSZ);

      auto concatParts = [&](const char *suffix, const char *header,
                             const char *tag, const char *outName) {
        std::string outFile = outputFile + outName;
        std::ofstream fout(outFile.c_str(),
                           std::ios::out | std::ios::trunc | std::ios::binary);
        fout << header;
        size_t rows = 0;
        for (int p = 0; p < world.size(); p++) {
          std::ostringstream oss;
          oss << outputFile << suffix << p << ".part";
          std::string part = oss.str();
          std::ifstream fin(part.c_str(), std::ios::in | std::ios::binary);
          if (!fin) continue;
          while (fin) {
            fin.read(&buf[0], BUFSZ);
            std::streamsize got = fin.gcount();
            if (got > 0) {
              fout.write(&buf[0], got);
              for (std::streamsize k = 0; k < got; k++) {
                if (buf[k] == '\n') rows++;
              }
            }
          }
          fin.close();
          std::remove(part.c_str());
        }
        fout.close();
        printf("[%s] wrote %zu rows to %s\n", tag, rows, outFile.c_str());
      };

      if (outputFormat == "csv_flat") {
        concatParts("_flat.",
                    "event_idx,contingency,from_bus,to_bus,circuit_id,"
                    "p_from_mw,q_from_mvar,mva_from,rate_mva,loading_percent,"
                    "viol,v_from_pu,v_to_pu,ang_from_deg,ang_to_deg\n",
                    "csv_flat",
                    "_flat.csv");
      }
      if (outputFormat == "csv_delta") {
        concatParts("_delta.",
                    "event_idx,contingency,type,from_bus,to_bus,ckt,"
                    "base_kv_from,base_kv_to,area_from,area_to,base_rate_mva,cont_rate_mva,"
                    "base_p_mw,cont_p_mw,base_q_mvar,cont_q_mvar,"
                    "base_mva,cont_mva,base_loading_pct,cont_loading_pct,"
                    "v_from_base,v_from_cont,v_to_base,v_to_cont,"
                    "ang_from_base,ang_from_cont,ang_to_base,ang_to_cont,"
                    // ",cont_event_facility" here if re-enabling the column
                    "d_v_base,d_v_cont,d_angle_base,d_angle_cont\n",
                    "csv_delta",
                    "_delta.csv");
      }

      // Bus metadata sidecar (deduped by bus_id, first writer wins).
      std::string busFile = outputFile + "_buses.csv";
      std::ofstream bout(busFile.c_str(),
                         std::ios::out | std::ios::trunc | std::ios::binary);
      bout << "bus_id,bus_name,base_kv,area,zone,owner,"
              "area_name,zone_name,owner_name\n";
      std::set<int> seen_bus;
      size_t bus_rows = 0;
      for (int p = 0; p < world.size(); p++) {
        std::ostringstream oss;
        oss << outputFile << "_buses." << p << ".part";
        std::string part = oss.str();
        std::ifstream fin(part.c_str());
        if (!fin) continue;
        std::string line;
        while (std::getline(fin, line)) {
          if (line.empty()) continue;
          size_t comma = line.find(',');
          if (comma == std::string::npos) continue;
          int bus_id = std::atoi(line.substr(0, comma).c_str());
          if (seen_bus.insert(bus_id).second) {
            bout << line << "\n";
            bus_rows++;
          }
        }
        fin.close();
        std::remove(part.c_str());
      }
      bout.close();
      printf("[buses] wrote %zu rows to %s\n", bus_rows, busFile.c_str());
    }
  }
  // Concat per-rank _violations.<rank>.part into <outputFile>_violations.csv.
  // Runs for every outputFormat; ranks that emitted zero rows simply have no
  // .part file to include.
  if (violPart.is_open()) violPart.close();
  world.sync();
  if (world.rank() == 0) {
    const size_t BUFSZ = 1 << 20;
    std::vector<char> buf(BUFSZ);
    std::string outFile = outputFile + "_violations.csv";
    std::ofstream fout(outFile.c_str(),
                       std::ios::out | std::ios::trunc | std::ios::binary);
    fout << "event_idx,contingency,type,element,mva_or_vpu,rate_or_limit,"
            "loading_percent,base_mva,delta,severity\n";
    size_t rows = 0;
    for (int p = 0; p < world.size(); p++) {
      std::ostringstream oss;
      oss << outputFile << "_violations." << p << ".part";
      std::string part = oss.str();
      std::ifstream fin(part.c_str(), std::ios::in | std::ios::binary);
      if (!fin) continue;
      while (fin) {
        fin.read(&buf[0], BUFSZ);
        std::streamsize got = fin.gcount();
        if (got > 0) {
          fout.write(&buf[0], got);
          for (std::streamsize k = 0; k < got; k++) {
            if (buf[k] == '\n') rows++;
          }
        }
      }
      fin.close();
      std::remove(part.c_str());
    }
    fout.close();
    printf("[violations] wrote %zu rows to %s\n", rows, outFile.c_str());
  }

  // Aggregate skip count across ranks for diagnostics.
  if (outputFormat == "csv_delta") {
    long localSkip = static_cast<long>(deltaSkipCount);
    long totalSkip = localSkip;
    world.sum(&totalSkip, 1);
    if (world.rank() == 0 && totalSkip > 0) {
      printf("[csv_delta] %ld branch rows had no base-cache match\n",
             totalSkip);
    }
  }

  // Aggregate per-run summary across all ranks. Emit <outputFile>_summary.json
  // on rank 0. Counters and worst-of values follow the same shape commercial
  // tools use in their contingency reports.
  {
    // 0: total_ct  1: converged  2: cts_with_branch_viol  3: cts_with_voltage_viol
    // 4: islanded  5: no_slack   6: diverged                7: slack_overload
    long localCounters[8] = { 0 };
    if (!localContingencies.empty()) {
      // json/csv paths: authoritative per-ct list.
      for (size_t ci = 0; ci < localContingencies.size(); ci++) {
        const gridpack::utility::ContingencyResult &ct = localContingencies[ci];
        localCounters[0] += 1;
        if (ct.solution.convergence.converged) localCounters[1] += 1;
      }
    } else {
      // text/csv_flat/csv_delta: count from convergence rows.
      localCounters[0] = static_cast<long>(localConvRows.size());
      long conv = 0;
      for (size_t i = 0; i < localConvRows.size(); i++) {
        if (localConvRows[i].status == "OK") conv++;
      }
      localCounters[1] = conv;
    }
    localCounters[2] = static_cast<long>(ctsWithBranchViol.size());
    localCounters[3] = static_cast<long>(ctsWithVoltageViol.size());
    // Per-status breakdown from convergence rows (populated in every mode).
    for (size_t i = 0; i < localConvRows.size(); i++) {
      const std::string &st = localConvRows[i].status;
      if      (st == "ISLANDED")       localCounters[4] += 1;
      else if (st == "NO_SLACK")       localCounters[5] += 1;
      else if (st == "DIVERGED")       localCounters[6] += 1;
      else if (st == "SLACK_OVERLOAD") localCounters[7] += 1;
    }
    long totalCounters[8] = { 0 };
    for (int i = 0; i < 8; i++) totalCounters[i] = localCounters[i];
    world.sum(&totalCounters[0], 8);
    long localViolRows  = static_cast<long>(violRowCount);
    long totalViolRows  = localViolRows;
    world.sum(&totalViolRows, 1);

    struct WorstBranchWire {
      double loading_pct;
      int from, to;
      char ckt[4];
      char ct_name[32];
    };
    struct WorstVoltageWire {
      double v_pu;
      double dev_pu;
      int bus_id;
      char ct_name[32];
    };
    WorstBranchWire wbLo;
    wbLo.loading_pct = worstBr.loading_pct;
    wbLo.from = worstBr.from; wbLo.to = worstBr.to;
    std::strncpy(wbLo.ckt, worstBr.ckt.c_str(), 3); wbLo.ckt[3] = '\0';
    std::strncpy(wbLo.ct_name, worstBr.ct_name.c_str(), 31); wbLo.ct_name[31] = '\0';
    WorstVoltageWire wvLo;
    wvLo.v_pu = worstVLo.v_pu; wvLo.dev_pu = worstVLo.dev_pu;
    wvLo.bus_id = worstVLo.bus_id;
    std::strncpy(wvLo.ct_name, worstVLo.ct_name.c_str(), 31); wvLo.ct_name[31] = '\0';
    WorstVoltageWire wvHi;
    wvHi.v_pu = worstVHi.v_pu; wvHi.dev_pu = worstVHi.dev_pu;
    wvHi.bus_id = worstVHi.bus_id;
    std::strncpy(wvHi.ct_name, worstVHi.ct_name.c_str(), 31); wvHi.ct_name[31] = '\0';
    MPI_Comm mpi_comm = static_cast<MPI_Comm>(world);
    if (world.rank() == 0) {
      for (int p = 1; p < world.size(); p++) {
        WorstBranchWire  otherB;
        WorstVoltageWire otherLo, otherHi;
        MPI_Recv(&otherB,  sizeof(otherB),  MPI_BYTE, p, 30, mpi_comm, MPI_STATUS_IGNORE);
        MPI_Recv(&otherLo, sizeof(otherLo), MPI_BYTE, p, 31, mpi_comm, MPI_STATUS_IGNORE);
        MPI_Recv(&otherHi, sizeof(otherHi), MPI_BYTE, p, 32, mpi_comm, MPI_STATUS_IGNORE);
        if (otherB.loading_pct > wbLo.loading_pct) wbLo = otherB;
        if (otherLo.dev_pu < wvLo.dev_pu)          wvLo = otherLo;
        if (otherHi.dev_pu > wvHi.dev_pu)          wvHi = otherHi;
      }
    } else {
      MPI_Send(&wbLo, sizeof(wbLo), MPI_BYTE, 0, 30, mpi_comm);
      MPI_Send(&wvLo, sizeof(wvLo), MPI_BYTE, 0, 31, mpi_comm);
      MPI_Send(&wvHi, sizeof(wvHi), MPI_BYTE, 0, 32, mpi_comm);
    }

    // Gather per-ct PI/VPI/Vdev/rosters + per-ct worst-single-element data.
    // Blob line: name pi vpi vdev worstLoad worstVdevAbs worstVpu br_flag v_flag
    std::map<std::string, double> aggPi;
    std::map<std::string, double> aggVpi;
    std::map<std::string, double> aggVdev;
    std::map<std::string, double> aggWorstLoading;
    std::map<std::string, double> aggWorstVdev;    // unsigned max |v-limit|
    std::map<std::string, double> aggWorstVpu;     // v_pu that produced aggWorstVdev
    std::set<std::string> aggBranchViol;
    std::set<std::string> aggVoltageViol;
    {
      std::ostringstream localOut;
      localOut << std::setprecision(10);
      std::set<std::string> names;
      for (std::map<std::string,double>::const_iterator it = ctPi.begin();
           it != ctPi.end(); ++it) names.insert(it->first);
      for (std::map<std::string,double>::const_iterator it = ctVpi.begin();
           it != ctVpi.end(); ++it) names.insert(it->first);
      for (std::map<std::string,double>::const_iterator it = ctVdev.begin();
           it != ctVdev.end(); ++it) names.insert(it->first);
      for (std::set<std::string>::const_iterator it = ctsWithBranchViol.begin();
           it != ctsWithBranchViol.end(); ++it) names.insert(*it);
      for (std::set<std::string>::const_iterator it = ctsWithVoltageViol.begin();
           it != ctsWithVoltageViol.end(); ++it) names.insert(*it);
      auto lookup = [](const std::map<std::string,double> &m,
                       const std::string &k) -> double {
        std::map<std::string,double>::const_iterator it = m.find(k);
        return (it == m.end()) ? 0.0 : it->second;
      };
      for (std::set<std::string>::const_iterator it = names.begin();
           it != names.end(); ++it) {
        const std::string &nm = *it;
        double pi    = lookup(ctPi,  nm);
        double vpi   = lookup(ctVpi, nm);
        double vdev  = lookup(ctVdev, nm);
        double wL    = lookup(ctWorstLoading, nm);
        double wDabs = lookup(ctWorstVdev, nm);
        double wV    = lookup(ctWorstVpu,  nm);
        int brFlag = ctsWithBranchViol.count(nm)  ? 1 : 0;
        int vFlag  = ctsWithVoltageViol.count(nm) ? 1 : 0;
        localOut << nm << '\t' << pi << '\t' << vpi << '\t' << vdev << '\t'
                 << wL << '\t' << wDabs << '\t' << wV << '\t'
                 << brFlag << '\t' << vFlag << '\n';
      }
      std::string localBlob = localOut.str();
      auto splitTabs = [](const std::string &line,
                          std::vector<std::string> &out) {
        out.clear();
        size_t pos = 0;
        while (pos <= line.size()) {
          size_t t = line.find('\t', pos);
          if (t == std::string::npos) {
            out.push_back(line.substr(pos));
            break;
          }
          out.push_back(line.substr(pos, t - pos));
          pos = t + 1;
        }
      };
      auto absorb = [&](const std::string &blob) {
        size_t pos = 0;
        std::vector<std::string> fields;
        while (pos < blob.size()) {
          size_t eol = blob.find('\n', pos);
          if (eol == std::string::npos) break;
          std::string line = blob.substr(pos, eol - pos);
          pos = eol + 1;
          splitTabs(line, fields);
          if (fields.size() < 9) continue;
          const std::string &nm = fields[0];
          double pi    = std::atof(fields[1].c_str());
          double vpi   = std::atof(fields[2].c_str());
          double vdev  = std::atof(fields[3].c_str());
          double wL    = std::atof(fields[4].c_str());
          double wDabs = std::atof(fields[5].c_str());
          double wV    = std::atof(fields[6].c_str());
          int brF      = std::atoi(fields[7].c_str());
          int vF       = std::atoi(fields[8].c_str());
          if (pi   != 0.0) aggPi[nm]   += pi;
          if (vpi  != 0.0) aggVpi[nm]  += vpi;
          if (vdev != 0.0) aggVdev[nm] += vdev;
          if (wL > aggWorstLoading[nm])                aggWorstLoading[nm] = wL;
          if (wDabs > aggWorstVdev[nm]) { aggWorstVdev[nm] = wDabs; aggWorstVpu[nm] = wV; }
          if (brF) aggBranchViol.insert(nm);
          if (vF)  aggVoltageViol.insert(nm);
        }
      };
      if (world.rank() == 0) {
        absorb(localBlob);
        for (int p = 1; p < world.size(); p++) {
          int len = 0;
          MPI_Recv(&len, 1, MPI_INT, p, 33, mpi_comm, MPI_STATUS_IGNORE);
          std::string remote;
          remote.resize(len);
          if (len > 0) {
            MPI_Recv(&remote[0], len, MPI_CHAR, p, 34, mpi_comm,
                     MPI_STATUS_IGNORE);
          }
          absorb(remote);
        }
      } else {
        int len = static_cast<int>(localBlob.size());
        MPI_Send(&len, 1, MPI_INT, 0, 33, mpi_comm);
        if (len > 0) {
          MPI_Send(const_cast<char*>(localBlob.c_str()), len, MPI_CHAR, 0, 34,
                   mpi_comm);
        }
      }
    }

    if (world.rank() == 0) {
      std::string sumFile = outputFile + "_summary.json";
      std::ofstream sout(sumFile.c_str());
      sout << std::fixed;
      sout << "{\n";
      sout << "  \"total_contingencies\": "     << totalCounters[0] << ",\n";
      sout << "  \"converged\": "               << totalCounters[1] << ",\n";
      sout << "  \"diverged\": "                << (totalCounters[0] - totalCounters[1]) << ",\n";
      // Per-status breakdown of the diverged bucket. Sum equals `diverged`.
      sout << "  \"islanded\": "                << totalCounters[4] << ",\n";
      sout << "  \"no_slack\": "                << totalCounters[5] << ",\n";
      sout << "  \"solver_diverged\": "         << totalCounters[6] << ",\n";
      sout << "  \"slack_overload\": "          << totalCounters[7] << ",\n";
      sout << "  \"with_branch_violation\": "   << totalCounters[2] << ",\n";
      sout << "  \"with_voltage_violation\": "  << totalCounters[3] << ",\n";
      sout << "  \"worst_loading\": ";
      if (wbLo.loading_pct > 0.0) {
        sout << "{\"contingency\": \"" << wbLo.ct_name << "\""
             << ", \"from_bus\": " << wbLo.from
             << ", \"to_bus\": "   << wbLo.to
             << ", \"circuit_id\": \"" << wbLo.ckt << "\""
             << ", \"loading_percent\": " << std::setprecision(2) << wbLo.loading_pct
             << "},\n";
      } else {
        sout << "null,\n";
      }
      sout << "  \"worst_voltage_low\": ";
      if (wvLo.dev_pu < 0.0) {
        sout << "{\"contingency\": \"" << wvLo.ct_name << "\""
             << ", \"bus_id\": " << wvLo.bus_id
             << ", \"v_pu\": "          << std::setprecision(6) << wvLo.v_pu
             << ", \"deviation_pu\": "  << std::setprecision(6) << wvLo.dev_pu
             << "},\n";
      } else {
        sout << "null,\n";
      }
      sout << "  \"worst_voltage_high\": ";
      if (wvHi.dev_pu > 0.0) {
        sout << "{\"contingency\": \"" << wvHi.ct_name << "\""
             << ", \"bus_id\": " << wvHi.bus_id
             << ", \"v_pu\": "         << std::setprecision(6) << wvHi.v_pu
             << ", \"deviation_pu\": " << std::setprecision(6) << wvHi.dev_pu
             << "},\n";
      } else {
        sout << "null,\n";
      }
      auto emitNameArray = [&](const char *field,
                               const std::set<std::string> &names) {
        sout << "  \"" << field << "\": [";
        int emitted = 0;
        for (std::set<std::string>::const_iterator it = names.begin();
             it != names.end(); ++it, ++emitted) {
          if (emitted) sout << ", ";
          sout << "\"" << *it << "\"";
        }
        sout << "],\n";
      };
      emitNameArray("contingencies_with_branch_violation", aggBranchViol);
      emitNameArray("contingencies_with_voltage_violation", aggVoltageViol);
      // top_severe_contingencies: Group A (any violation) first, sorted by
      // worst-single-element severity; Group B (no violations) after,
      // sorted by composite_pi. Combined list capped at topN.
      {
        std::set<std::string> ctSet;
        for (std::map<std::string,double>::const_iterator it = aggPi.begin();
             it != aggPi.end(); ++it) ctSet.insert(it->first);
        for (std::map<std::string,double>::const_iterator it = aggVpi.begin();
             it != aggVpi.end(); ++it) ctSet.insert(it->first);
        for (std::set<std::string>::const_iterator it = aggBranchViol.begin();
             it != aggBranchViol.end(); ++it) ctSet.insert(*it);
        for (std::set<std::string>::const_iterator it = aggVoltageViol.begin();
             it != aggVoltageViol.end(); ++it) ctSet.insert(*it);
        auto agg_get = [](const std::map<std::string,double> &m,
                          const std::string &k) -> double {
          std::map<std::string,double>::const_iterator it = m.find(k);
          return (it == m.end()) ? 0.0 : it->second;
        };
        struct SevRow {
          std::string name;
          bool violated;
          double sortKey;   // Group A: worst-single severity; Group B: composite_pi
          double composite, branchPi, voltagePi;
          double worstLoading;
          double worstVpu, worstVdev;   // worstVdev is unsigned |v - limit|
        };
        std::vector<SevRow> groupA, groupB;
        for (std::set<std::string>::const_iterator it = ctSet.begin();
             it != ctSet.end(); ++it) {
          SevRow r;
          r.name = *it;
          r.branchPi   = agg_get(aggPi, *it);
          r.voltagePi  = agg_get(aggVpi, *it);
          r.composite  = piBranchWeight * r.branchPi + piVoltageWeight * r.voltagePi;
          r.worstLoading = agg_get(aggWorstLoading, *it);
          r.worstVdev  = agg_get(aggWorstVdev, *it);
          r.worstVpu   = agg_get(aggWorstVpu,  *it);
          bool hasBr = aggBranchViol.count(*it) > 0;
          bool hasV  = aggVoltageViol.count(*it) > 0;
          r.violated = hasBr || hasV;
          if (r.violated) {
            double s = 0.0;
            if (r.worstLoading > 100.0) s = std::max(s, r.worstLoading - 100.0);
            if (r.worstVdev    > 0.0)   s = std::max(s, r.worstVdev * 1000.0);
            r.sortKey = s;
            groupA.push_back(r);
          } else {
            r.sortKey = r.composite;
            if (r.sortKey > 0.0) groupB.push_back(r);
          }
        }
        auto sevCmp = [](const SevRow &a, const SevRow &b) {
          if (a.sortKey != b.sortKey) return a.sortKey > b.sortKey;
          return a.name < b.name;
        };
        std::sort(groupA.begin(), groupA.end(), sevCmp);
        std::sort(groupB.begin(), groupB.end(), sevCmp);
        auto emitRow = [&](const SevRow &r, bool first) {
          if (!first) sout << ",\n    ";
          else        sout << "\n    ";
          bool hasBr = aggBranchViol.count(r.name) > 0;
          bool hasV  = r.worstVdev > 0.0;
          sout << "{\"contingency\": \"" << r.name << "\""
               << ", \"composite_pi\": " << std::setprecision(6) << r.composite
               << ", \"worst_branch_loading_percent\": ";
          if (hasBr) sout << std::setprecision(2) << r.worstLoading;
          else       sout << "null";
          sout << ", \"worst_voltage_pu\": ";
          if (hasV) sout << std::setprecision(6) << r.worstVpu;
          else      sout << "null";
          sout << ", \"worst_voltage_deviation_pu\": ";
          if (hasV) sout << std::setprecision(6) << r.worstVdev;
          else      sout << "null";
          sout << ", \"has_branch_violation\": "
               << (hasBr ? "true" : "false")
               << ", \"has_voltage_violation\": "
               << (aggVoltageViol.count(r.name) ? "true" : "false")
               << "}";
        };
        sout << "  \"top_severe_contingencies\": [";
        int emitted = 0;
        for (size_t i = 0; i < groupA.size() && emitted < topN; ++i, ++emitted) {
          emitRow(groupA[i], emitted == 0);
        }
        for (size_t i = 0; i < groupB.size() && emitted < topN; ++i, ++emitted) {
          emitRow(groupB[i], emitted == 0);
        }
        sout << (emitted ? "\n  ]\n" : "]\n");
      }
      sout << "}\n";
      sout.close();
      printf("[summary] wrote %s (%ld contingencies, %ld converged, "
             "%ld with branch violations, %ld with voltage violations)\n",
             sumFile.c_str(),
             totalCounters[0], totalCounters[1],
             totalCounters[2], totalCounters[3]);
    }
  }

  // Print statistics from task manager describing the number of tasks performed
  // per processor
  taskmgr.printStats();

  // Sync GA before MPI collectives to flush any pending one-sided operations
  world.sync();

  // Export CA results to JSON or CSV.
  // Only rank 0 writes output files. Non-zero ranks send their serialized
  // data to rank 0 using point-to-point MPI send/recv.
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

    // Gather all JSON fragments on rank 0 using point-to-point send/recv
    MPI_Comm mpi_comm = static_cast<MPI_Comm>(world);
    std::vector<std::string> allFragments(world.size());
    allFragments[0] = localJson;  // rank 0's own data
    if (world.rank() == 0) {
      for (int p = 1; p < world.size(); p++) {
        int len;
        MPI_Recv(&len, 1, MPI_INT, p, 0, mpi_comm, MPI_STATUS_IGNORE);
        allFragments[p].resize(len);
        if (len > 0) {
          MPI_Recv(&allFragments[p][0], len, MPI_CHAR, p, 1, mpi_comm,
                   MPI_STATUS_IGNORE);
        }
      }
    } else {
      int len = static_cast<int>(localJson.size());
      MPI_Send(&len, 1, MPI_INT, 0, 0, mpi_comm);
      if (len > 0) {
        MPI_Send(const_cast<char*>(localJson.c_str()), len, MPI_CHAR, 0, 1,
                 mpi_comm);
      }
    }

    // Rank 0 writes the final JSON file
    if (world.rank() == 0) {
      std::string jsonFile = outputFile + ".json";
      std::ofstream jout(jsonFile.c_str());
      gridpack::utility::ResultsExporter::writeCAJSONHeader(jout,
          baseCaseResults);
      bool firstFragment = true;
      for (int p = 0; p < world.size(); p++) {
        if (!allFragments[p].empty()) {
          if (!firstFragment) jout << ",\n";
          jout << allFragments[p];
          firstFragment = false;
        }
      }
      jout << "\n";
      gridpack::utility::ResultsExporter::writeCAJSONFooter(jout);
      jout.close();
    }
  }

  if (outputFormat == "csv") {
    // Convergence is emitted by the universal sidecar block below.
    std::ostringstream localBus, localBranch, localGen;
    localBus << std::fixed;
    localBranch << std::fixed;
    localGen << std::fixed;
    // Column layout must match the headers written by
    // ResultsExporter::writePFCSV for the base case.
    for (size_t ci = 0; ci < localContingencies.size(); ci++) {
      const gridpack::utility::ContingencyResult& ct = localContingencies[ci];
      const gridpack::utility::PowerFlowResults& r = ct.solution;
      for (size_t bi = 0; bi < r.buses.size(); bi++) {
        const gridpack::utility::BusResult& b = r.buses[bi];
        localBus << ct.name << ","
           << b.busId << "," << b.type << ","
           << b.area << "," << b.zone << ","
           << std::setprecision(2) << b.baseKV << ","
           << std::setprecision(6) << b.voltage << ","
           << std::setprecision(6) << b.angle << ","
           << std::setprecision(4) << b.pInjection << ","
           << std::setprecision(4) << b.qInjection << ","
           << std::setprecision(4) << b.pLoad << ","
           << std::setprecision(4) << b.qLoad << ","
           << std::setprecision(4) << b.pGen << ","
           << std::setprecision(4) << b.qGen << ","
           << std::setprecision(4) << b.shuntMvar << "\n";
      }
      for (size_t bi = 0; bi < r.branches.size(); bi++) {
        const gridpack::utility::BranchResult& br = r.branches[bi];
        localBranch << ct.name << ","
           << br.fromBus << "," << br.toBus << ","
           << br.circuitId << ","
           << std::setprecision(4) << br.pFrom << ","
           << std::setprecision(4) << br.qFrom << ","
           << std::setprecision(4) << br.pTo << ","
           << std::setprecision(4) << br.qTo << ","
           << std::setprecision(4) << br.pLoss << ","
           << std::setprecision(4) << br.qLoss << ","
           << std::setprecision(4) << br.mvaFrom << ","
           << std::setprecision(4) << br.mvaTo << ","
           << std::setprecision(4) << br.rateA << ","
           << std::setprecision(4) << br.rateSelected << ","
           << std::setprecision(2) << br.loadingPercent << "\n";
      }
      for (size_t gi = 0; gi < r.generators.size(); gi++) {
        const gridpack::utility::GeneratorResult& g = r.generators[gi];
        localGen << ct.name << ","
           << g.busId << "," << g.genId << ","
           << std::setprecision(4) << g.pGen << ","
           << std::setprecision(4) << g.qGen << ","
           << std::setprecision(4) << g.qMax << ","
           << std::setprecision(4) << g.qMin << ","
           << std::setprecision(6) << g.voltageSetpoint << ","
           << g.status << "\n";
      }
    }

    // Gather all CSV fragments on rank 0 using point-to-point send/recv
    MPI_Comm mpi_comm = static_cast<MPI_Comm>(world);
    std::vector<std::string> allBus(world.size()), allBranch(world.size());
    std::vector<std::string> allGen(world.size());
    allBus[0] = localBus.str();
    allBranch[0] = localBranch.str();
    allGen[0] = localGen.str();
    if (world.rank() == 0) {
      for (int p = 1; p < world.size(); p++) {
        int lens[3];
        MPI_Recv(lens, 3, MPI_INT, p, 0, mpi_comm, MPI_STATUS_IGNORE);
        allBus[p].resize(lens[0]);
        allBranch[p].resize(lens[1]);
        allGen[p].resize(lens[2]);
        if (lens[0] > 0)
          MPI_Recv(&allBus[p][0], lens[0], MPI_CHAR, p, 1, mpi_comm,
                   MPI_STATUS_IGNORE);
        if (lens[1] > 0)
          MPI_Recv(&allBranch[p][0], lens[1], MPI_CHAR, p, 2, mpi_comm,
                   MPI_STATUS_IGNORE);
        if (lens[2] > 0)
          MPI_Recv(&allGen[p][0], lens[2], MPI_CHAR, p, 3, mpi_comm,
                   MPI_STATUS_IGNORE);
      }
    } else {
      std::string sBus = localBus.str(), sBranch = localBranch.str();
      std::string sGen = localGen.str();
      int lens[3] = {(int)sBus.size(), (int)sBranch.size(), (int)sGen.size()};
      MPI_Send(lens, 3, MPI_INT, 0, 0, mpi_comm);
      if (lens[0] > 0)
        MPI_Send(const_cast<char*>(sBus.c_str()), lens[0], MPI_CHAR, 0, 1,
                 mpi_comm);
      if (lens[1] > 0)
        MPI_Send(const_cast<char*>(sBranch.c_str()), lens[1], MPI_CHAR, 0, 2,
                 mpi_comm);
      if (lens[2] > 0)
        MPI_Send(const_cast<char*>(sGen.c_str()), lens[2], MPI_CHAR, 0, 3,
                 mpi_comm);
    }

    // Rank 0 writes the CSV files
    if (world.rank() == 0) {
      // Write base case first (creates files with headers)
      gridpack::utility::ResultsExporter::writePFCSV(outputFile,
          baseCaseResults, "base_case");
      // Append contingency data from all processes
      {
        std::ofstream out((outputFile + "_buses.csv").c_str(), std::ios::app);
        for (size_t p = 0; p < allBus.size(); p++) out << allBus[p];
      }
      {
        std::ofstream out((outputFile + "_branches.csv").c_str(), std::ios::app);
        for (size_t p = 0; p < allBranch.size(); p++) out << allBranch[p];
      }
      {
        std::ofstream out((outputFile + "_generators.csv").c_str(), std::ios::app);
        for (size_t p = 0; p < allGen.size(); p++) out << allGen[p];
      }
    }
  }

  // Universal convergence sidecar: gather, sort by event_idx, write.
  if (emitConv) {
    auto formatRow = [](std::ostringstream &os, const ConvRow &r) {
      // Derive converged from status so ISLANDED/NO_SLACK cases -- where solve()
      // was never entered and pf_app.getConvergence() returns the previous
      // case's stale value -- read as false, matching _summary.json's
      // diverged=total-converged accounting.
      os << r.event_idx << ","
         << r.name << ","
         << r.type << ","
         << ((r.status == "OK") ? "true" : "false") << ","
         << r.cs.iterations << ","
         << std::scientific << r.cs.finalTolerance << ","
         << std::fixed
         << r.cs.finalMismatch.maxPBus << ","
         << std::setprecision(4) << r.cs.finalMismatch.maxPMismatch << ","
         << r.cs.finalMismatch.maxQBus << ","
         << std::setprecision(4) << r.cs.finalMismatch.maxQMismatch << ","
         << r.status << "\n";
    };

    std::vector<int> idx;
    std::ostringstream localStream;
    localStream << std::fixed;
    std::vector<int> localOffsets;
    localOffsets.reserve(localConvRows.size() + 1);
    for (size_t i = 0; i < localConvRows.size(); i++) {
      localOffsets.push_back(static_cast<int>(localStream.tellp()));
      formatRow(localStream, localConvRows[i]);
      idx.push_back(localConvRows[i].event_idx);
    }
    localOffsets.push_back(static_cast<int>(localStream.tellp()));
    std::string localStr = localStream.str();

    MPI_Comm conv_comm = static_cast<MPI_Comm>(world);
    if (world.rank() == 0) {
      std::vector<std::pair<int, std::string> > all;
      for (size_t i = 0; i < idx.size(); i++) {
        std::string row = localStr.substr(localOffsets[i],
                                          localOffsets[i+1] - localOffsets[i]);
        all.push_back(std::make_pair(idx[i], row));
      }
      for (int p = 1; p < world.size(); p++) {
        int n = 0;
        MPI_Recv(&n, 1, MPI_INT, p, 10, conv_comm, MPI_STATUS_IGNORE);
        if (n <= 0) continue;
        std::vector<int> remIdx(n), remOff(n + 1);
        MPI_Recv(&remIdx[0], n, MPI_INT, p, 11, conv_comm, MPI_STATUS_IGNORE);
        MPI_Recv(&remOff[0], n + 1, MPI_INT, p, 12, conv_comm,
                 MPI_STATUS_IGNORE);
        int total = remOff[n];
        std::string buf(total, '\0');
        if (total > 0) {
          MPI_Recv(&buf[0], total, MPI_CHAR, p, 13, conv_comm,
                   MPI_STATUS_IGNORE);
        }
        for (int i = 0; i < n; i++) {
          all.push_back(std::make_pair(
              remIdx[i],
              buf.substr(remOff[i], remOff[i+1] - remOff[i])));
        }
      }
      std::sort(all.begin(), all.end());
      std::string convFile = outputFile + "_convergence.csv";
      std::ofstream cout(convFile.c_str(),
                         std::ios::out | std::ios::trunc);
      cout << "event_idx,contingency,type,converged,iterations,"
              "final_tolerance,max_p_bus,max_p_mismatch,max_q_bus,"
              "max_q_mismatch,status_code\n";
      for (size_t i = 0; i < all.size(); i++) cout << all[i].second;
      cout.close();
      printf("[convergence] wrote %zu rows to %s\n",
             all.size(), convFile.c_str());
    } else {
      int n = static_cast<int>(idx.size());
      MPI_Send(&n, 1, MPI_INT, 0, 10, conv_comm);
      if (n > 0) {
        MPI_Send(&idx[0], n, MPI_INT, 0, 11, conv_comm);
        MPI_Send(&localOffsets[0], n + 1, MPI_INT, 0, 12, conv_comm);
        int total = localOffsets[n];
        if (total > 0) {
          MPI_Send(const_cast<char*>(localStr.c_str()), total, MPI_CHAR, 0, 13,
                   conv_comm);
        }
      }
    }
  }

  // Print out statistics on contingencies
  if (write_stats) {
    int t_stats = timer->createCategory("Write Statistics");
    timer->start(t_stats);
    vmag_stats->writeMeanAndRMS("vmag.txt",1,false);
    vmag_stats->writeMinAndMax("vmag_mm.txt",1,false);
    if (check_Qlim) vmag_stats->writeMaskValueCount("pq_change_cnt.txt",2,false);
    vang_stats->writeMeanAndRMS("vang.txt",1,false);
    vang_stats->writeMinAndMax("vang_mm.txt",1,false);
    pgen_stats->writeMeanAndRMS("pgen.txt",1);
    pgen_stats->writeMinAndMax("pgen_mm.txt",1);
    qgen_stats->writeMeanAndRMS("qgen.txt",1);
    qgen_stats->writeMinAndMax("qgen_mm.txt",1);
    pflow_stats->writeMeanAndRMS("pflow.txt",1);
    pflow_stats->writeMinAndMax("pflow_mm.txt",1);
    pflow_stats->writeMaskValueCount("line_flt_cnt.txt",2);
    qflow_stats->writeMeanAndRMS("qflow.txt",1);
    qflow_stats->writeMinAndMax("qflow_mm.txt",1);
    perf_stats->writeMinAndMax("perf_mm.txt",1);
    perf_stats->sumColumnValues("perf_sum.txt",1);
    timer->stop(t_stats);
  }
  timer->stop(t_total);
  // If all processors executed at least one task, then print out timing
  // statistics (this printout does not work if some processors do not define
  // all timing variables)
  if (events.size()*grp_size >= world.size()) {
    timer->dump();
  }
}

