/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsf_main.cpp
 * @author Shuangshuang Jin
 * @date   2016-07-14 14:23:30 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/math/math.hpp"
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
    ///printf("Step0 bus%d mag = %f\n", i+1, rval);
    pfData->getValue("BUS_PF_VANG",&rval);
    dsData->setValue(BUS_VOLTAGE_ANG,rval);
    int ngen = 0;
    if (pfData->getValue(GENERATOR_NUMBER, &ngen)) {
      int j;
      for (j=0; j<ngen; j++) {
        pfData->getValue("GENERATOR_PF_PGEN",&rval,j);
        dsData->setValue(GENERATOR_PG,rval,j);
        //printf("save PGEN: %f\n", rval);
        pfData->getValue("GENERATOR_PF_QGEN",&rval,j);
        dsData->setValue(GENERATOR_QG,rval,j);
        //printf("save QGEN: %f\n", rval);
      }
    }
  }
}

// Calling program for the dynamis simulation applications

int
main(int argc, char **argv)
{
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);

  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  // Intialize Math libraries
  gridpack::math::Initialize(&argc,&argv);

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

    // setup and run powerflow calculation
    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config->getCursor("Configuration.Powerflow");
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
   
    // setup and run dynamic simulation calculation
    boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork>
      ds_network(new gridpack::dynamic_simulation::DSFullNetwork(world));
    gridpack::dynamic_simulation::DSFullApp ds_app;
    pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
      gridpack::dynamic_simulation::DSFullBranch>(ds_network);

    // transfer results from PF calculation to DS calculation
    transferPFtoDS(pf_network, ds_network); 

    // read in faults from input file
    //gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config->getCursor("Configuration.Dynamic_simulation");
    std::vector<gridpack::dynamic_simulation::Event> faults;
    faults = ds_app.getFaults(cursor);

    // run dynamic simulation
    ds_app.setNetwork(ds_network, config);
    //ds_app.readNetwork(ds_network,config);
    ds_app.readGenerators();
    //printf("ds_app.initialize:\n");
    ds_app.initialize();
    ds_app.setGeneratorWatch();
    //printf("gen ID:	mac_ang_s0	mac_spd_s0	pmech	pelect\n");
    //printf("Step	time:	bus_id	mac_ang_s1	mac_spd_s1\n");
    //printf("ds_app.solve:\n");
    ds_app.solve(faults[0]);
    //ds_app.write();
    timer->stop(t_total);
    timer->dump();
  }

  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
}

