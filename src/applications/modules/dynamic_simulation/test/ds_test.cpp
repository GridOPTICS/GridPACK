/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_test.cpp
 * @author Bruce Palmer
 * @date   2016-07-14 14:23:29 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "gridpack/applications/modules/dynamic_simulation/ds_app_module.hpp"

/**
 * Transfer data from power flow to dynamic simulation
 * @param pf_network power flow network
 * @param ds_network dynamic simulation network
 */
void transferPFtoDS(
    boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network,
    boost::shared_ptr<gridpack::dynamic_simulation::DSNetwork>
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

// Test program for dynamic simulation module using reduced y-matrix

int
main(int argc, char **argv)
{
  gridpack::Environment env(argc,argv);

  if (1) {
    gridpack::parallel::Communicator world;

    // read configuration file
    gridpack::utility::Configuration *config =
      gridpack::utility::Configuration::configuration();
    if (argc >= 2 && argv[1] != NULL) {
      char inputfile[256];
      sprintf(inputfile,"%s",argv[1]);
      config->open(inputfile,world);
    } else {
      config->open("input.xml",world);
    }

    // create dynamic simulation network and app
    boost::shared_ptr<gridpack::dynamic_simulation::DSNetwork>
      ds_network(new gridpack::dynamic_simulation::DSNetwork(world));
    gridpack::dynamic_simulation::DSAppModule ds_app;

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
      pf_network->clone<gridpack::dynamic_simulation::DSBus,
        gridpack::dynamic_simulation::DSBranch>(ds_network);
      // transfer results from PF calculation to DS calculation
      transferPFtoDS(pf_network, ds_network);
      //             // set network for dynamic simulation
      ds_app.setNetwork(ds_network, config);
    } else {
      // read in network
      ds_app.readNetwork(ds_network,config);
    }

    ds_app.readGenerators();
    cursor = config->getCursor("Configuration.Dynamic_simulation");
    std::vector<gridpack::dynamic_simulation::DSBranch::Event> faults;
    faults = ds_app.getFaults(cursor);
    ds_app.initialize();
    ds_app.setGeneratorWatch();
    ds_app.solve(faults[0]);
    ds_app.write();
  }

  return 0;
}

