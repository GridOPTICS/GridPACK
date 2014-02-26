/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ca_driver.cpp
 * @author Bruce Palmer
 * @date   February 10, 2014
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#include "gridpack/parallel/task_manager.hpp"
#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/math/linear_solver.hpp"
#include "gridpack/math/linear_matrix_solver.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/applications/contingency_analysis/ca_driver.hpp"
#include "gridpack/applications/contingency_analysis/ca_app.hpp"

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
std::vector<gridpack::contingency_analysis::Contingency>
  gridpack::contingency_analysis::CADriver::getContingencies(
      gridpack::utility::Configuration::ChildCursors contingencies)
{
  std::vector<gridpack::contingency_analysis::Contingency> ret;
  int size = contingencies.size();
  int idx;
  for (idx = 0; idx < size; idx++) {
    std::string ca_type;
    contingencies[idx]->get("contingencyType",&ca_type);
    std::string ca_name;
    contingencies[idx]->get("contingencyName",&ca_name);
    if (ca_type == "Line") {
      std::string buses;
      contingencies[idx]->get("contingencyLineBuses",&buses);
      std::string names;
      contingencies[idx]->get("contingencyLineNames",&names);
      // Parse buses to get bus IDs in contingency
      int ntok1;
      int ntok2;
      ntok1 = buses.find_first_not_of(' ',0);
      ntok2 = buses.find(' ',ntok1);
      std::vector<int> bus_ids;
      while (ntok1 != std::string::npos) {
        if (ntok2 == std::string::npos) ntok2 = buses.length();
        if (ntok2<=ntok1) break;
        int bus = atoi(buses.substr(ntok1,ntok2-ntok1).c_str());
        bus_ids.push_back(bus);
        ntok1 = buses.find_first_not_of(' ',ntok2);
        ntok2 = buses.find(' ',ntok1);
      }
      // Parse names to get line names
      ntok1 = names.find_first_not_of(' ',0);
      ntok2 = names.find(' ',ntok1);
      std::vector<std::string> line_names;
      while (ntok1 != std::string::npos) {
        if (ntok2 == std::string::npos) ntok2 = names.length();
        if (ntok2<=ntok1) break;
        std::string name = names.substr(ntok1,ntok2-ntok1).c_str();
        // Pad single character names with a blank
        if (name.length() == 1) {
          std::string tmp = " ";
          tmp.append(name);
          name = tmp;
        }
        line_names.push_back(name);
        ntok1 = names.find_first_not_of(' ',ntok2);
        ntok2 = names.find(' ',ntok1);
      }
      // Check to make sure we found everything
      if (bus_ids.size() == 2*line_names.size()) {
        gridpack::contingency_analysis::Contingency contingency;
        contingency.p_name = ca_name;
        contingency.p_id = idx + 1;
        contingency.p_type = Branch;
        int i;
        for (i = 0; i < line_names.size(); i++) {
          contingency.p_from.push_back(bus_ids[2*i]);
          contingency.p_to.push_back(bus_ids[2*i+1]);
          contingency.p_ckt.push_back(line_names[i]);
        }
        ret.push_back(contingency);
      }
    } else if (ca_type == "Generator") {
      std::string buses;
      contingencies[idx]->get("contingencyBuses",&buses);
      std::string gens;
      contingencies[idx]->get("contingencyGenerators",&gens);
      // Parse buses to get bus IDs in contingency
      int ntok1;
      int ntok2;
      ntok1 = buses.find_first_not_of(' ',0);
      ntok2 = buses.find(' ',ntok1);
      std::vector<int> bus_ids;
      while (ntok1 != std::string::npos) {
        if (ntok2 == std::string::npos) ntok2 = buses.length();
        if (ntok2<=ntok1) break;
        int bus = atoi(buses.substr(ntok1,ntok2-ntok1).c_str());
        bus_ids.push_back(bus);
        ntok1 = buses.find_first_not_of(' ',ntok2);
        ntok2 = buses.find(' ',ntok1);
      }
      // Parse gens to get generator IDs in contingency
      ntok1 = gens.find_first_not_of(' ',0);
      ntok2 = gens.find(' ',ntok1);
      std::vector<std::string> gen_ids;
      while (ntok1 != std::string::npos) {
        if (ntok2 == std::string::npos) ntok2 = gens.length();
        if (ntok2<=ntok1) break;
        std::string gen = gens.substr(ntok1,ntok2-ntok1);
        if (gen.length() == 1) gen.insert(gen.begin(),' ');
        gen_ids.push_back(gen);
        ntok1 = gens.find_first_not_of(' ',ntok2);
        ntok2 = gens.find(' ',ntok1);
      }
      // Check to make sure we found everything
      if (bus_ids.size() == gen_ids.size()) {
        gridpack::contingency_analysis::Contingency contingency;
        contingency.p_name = ca_name;
        contingency.p_id = idx + 1;
        contingency.p_type = Generator;
        int i;
        for (i = 0; i < bus_ids.size(); i++) {
          contingency.p_busid.push_back(bus_ids[i]);
          contingency.p_genid.push_back(gen_ids[i]);
        }
        ret.push_back(contingency);
      }
    }
  }
  return ret;
}

/**
 * Execute application
 */
void gridpack::contingency_analysis::CADriver::execute(int argc, char** argv)
{
  gridpack::parallel::Communicator world;

  // read configuration file
  gridpack::utility::Configuration *config = gridpack::utility::Configuration::configuration();
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
  cursor = config->getCursor("Configuration.Contingency_analysis");
  int grp_size;
  if (!cursor->get("groupSize",&grp_size)) {
    grp_size = 1;
  }
  gridpack::parallel::Communicator task_comm = world.divide(grp_size);


  // get a list of contingencies
  cursor =
    config->getCursor("Configuration.Contingency_analysis.Contingencies");
  gridpack::utility::Configuration::ChildCursors contingencies;
  if (cursor) cursor->children(contingencies);
  std::vector<gridpack::contingency_analysis::Contingency>
    events = getContingencies(contingencies);
  if (world.rank() == 0) {
    int idx;
    for (idx = 0; idx < events.size(); idx++) {
      printf("Name: %s\n",events[idx].p_name.c_str());
      if (events[idx].p_type == Branch) {
        int nlines = events[idx].p_from.size();
        int j;
        for (j=0; j<nlines; j++) {
          printf(" Line: (to) %d (from) %d (line) \'%s\'\n",
              events[idx].p_from[j],events[idx].p_to[j],
              events[idx].p_ckt[j].c_str());
        }
      } else if (events[idx].p_type == Generator) {
        int nbus = events[idx].p_busid.size();
        int j;
        for (j=0; j<nbus; j++) {
          printf(" Generator: (bus) %d (generator ID) %s\n",
              events[idx].p_busid[j],events[idx].p_genid[j].c_str());
        }
      }
    }
  }

  // Create contingency applications on each task communicator
  gridpack::contingency_analysis::CAApp ca_app(task_comm);
  ca_app.init(argc,argv);

  // set up task manager
  gridpack::parallel::TaskManager taskmgr(world);
  int ntasks = contingencies.size();
  taskmgr.set(ntasks);

  // evaluate contingencies
  int task_id;
  while (taskmgr.nextTask(task_comm, &task_id)) {
    printf("Executing task %d\n",task_id);
    ca_app.execute(events[task_id]);
  }
}

