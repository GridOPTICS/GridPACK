/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_main.cpp
 * @author Shuangshuang Jin
 * @date   2016-07-14 14:30:38 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "gridpack/applications/modules/dynamic_simulation_full_y/dsf_app_module.hpp"

/**
 * Transfer data from power flow to dynamic simulation
 * @param pf_network power flow network
 * @param ds_network dynamic simulation network
 */
void transferPFtoDS(
    boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network,
    boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork>
    ds_network)
{
  int numBus = pf_network->numBuses();
  int i;
  gridpack::component::DataCollection *pfData;
  gridpack::component::DataCollection *dsData;
  double rval;
  for (i=0; i<numBus; i++) {
    pfData = pf_network->getBusData(i).get();
    dsData = ds_network->getBusData(i).get();
    pfData->getValue("BUS_PF_VMAG",&rval);
    dsData->setValue(BUS_VOLTAGE_MAG,rval);
    pfData->getValue("BUS_PF_VANG",&rval);
    dsData->setValue(BUS_VOLTAGE_ANG,rval);
    int ngen = 0;
    if (pfData->getValue(GENERATOR_NUMBER, &ngen)) {
      int j;
      for (j=0; j<ngen; j++) {
        pfData->getValue("GENERATOR_PF_PGEN",&rval,j);
        dsData->setValue(GENERATOR_PG,rval,j);
        pfData->getValue("GENERATOR_PF_QGEN",&rval,j);
        dsData->setValue(GENERATOR_QG,rval,j);
      }
    }
  }
}


// test program for dynamic simulation module using full y-matrix

int
main(int argc, char **argv)
{
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);

  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  // Intialize Math libraries
  gridpack::math::Initialize();

  if (1) {
    gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
    int t_total = timer->createCategory("Dynamic Simulation: Total Application");
    timer->start(t_total);

    gridpack::parallel::Communicator world;

    // read configuration file 
    int t_config = timer->createCategory("Dynamic Simulation: Config");
    timer->start(t_config);
    gridpack::utility::Configuration *config =
      gridpack::utility::Configuration::configuration();
    if (argc >= 2 && argv[1] != NULL) { 
      char inputfile[256]; 
      sprintf(inputfile,"%s",argv[1]);
      config->open(inputfile,world);
    } else {
      config->open("input.xml",world);
    }
    timer->stop(t_config);

    // create dynamic simulation network and app
    boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork>
      ds_network(new gridpack::dynamic_simulation::DSFullNetwork(world));
    gridpack::dynamic_simulation::DSFullApp ds_app;

    // setup optional powerflow calculation to initialize dynamic simulation
    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config->getCursor("Configuration.Powerflow");
    if (cursor) {
      bool useNonLinear = false;
      useNonLinear = cursor->get("UseNonLinear", useNonLinear);

      boost::shared_ptr<gridpack::powerflow::PFNetwork>
        pf_network(new gridpack::powerflow::PFNetwork(world));

      gridpack::powerflow::PFAppModule pf_app;
      pf_app.readNetwork(pf_network, config);
      pf_app.initialize();
      if (useNonLinear) {
        pf_app.nl_solve();
      } else {
        pf_app.solve();
      }
      pf_app.write();
      pf_app.saveData();
      pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
        gridpack::dynamic_simulation::DSFullBranch>(ds_network);
      // transfer results from PF calculation to DS calculation
      transferPFtoDS(pf_network, ds_network);
      // set network for dynamic simulation
      ds_app.setNetwork(ds_network, config);
    } else {
      // read in network
      cursor = config->getCursor("Configuration.Dynamic_simulation");
      std::string basefile;
      int filetype = gridpack::dynamic_simulation::PTI23;
      if (!cursor->get("networkConfiguration",&basefile)) {
        cursor->get("networkConfiguration_v33",&basefile);
        filetype = gridpack::dynamic_simulation::PTI33;
      }
      ds_app.readNetwork(ds_network,config,basefile.c_str(),filetype);
    }
    // read in information on generators
    ds_app.readGenerators();

    // read in faults from input file
    cursor = config->getCursor("Configuration.Dynamic_simulation");
    std::vector<gridpack::dynamic_simulation::DSFullBranch::Event> faults;
    faults = ds_app.getFaults(cursor);

    // run dynamic simulation
    ds_app.initialize();
    ds_app.setGeneratorWatch();
    ds_app.solve(faults[0]);
    // ds_app.write();
    timer->stop(t_total);
    timer->dump();
  }

  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
}

