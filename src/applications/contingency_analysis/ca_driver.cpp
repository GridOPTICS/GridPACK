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
      gridpack::utility::Configuration::CursorPtr &cursor)
{
  std::vector<gridpack::contingency_analysis::Contingency> ret;
  // TODO: add code here
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
  std::vector<gridpack::contingency_analysis::Contingency>
    contingencies = getContingencies(cursor);

  // set up task manager
  gridpack::parallel::TaskManager taskmgr(world);
  int ntasks = contingencies.size();
  taskmgr.set(ntasks);

  // evaluate contingencies
  int task_id;
  while (taskmgr.nextTask(task_comm, &task_id)) {
    // ca_app.execute(contingencies[task_id]);
  }
}

