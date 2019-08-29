/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_main.cpp
 * @author Shuangshuang Jin
 * @date   Feb 04, 2015
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/math/math.hpp"
#include "dsf_app_module.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"

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
  int ival;
  double temp;
  double pi = 4.0*atan(1.0);
  int itmp = 0;
  for (i=0; i<numBus; i++) {
    pfData = pf_network->getBusData(i).get();
    dsData = ds_network->getBusData(i).get();
    pfData->getValue("BUS_PF_VMAG",&rval);
    dsData->setValue(BUS_VOLTAGE_MAG,rval);
    temp = rval;
    pfData->getValue("BUS_PF_VANG",&rval);
    dsData->setValue(BUS_VOLTAGE_ANG,rval);
    rval = rval * pi/180.0;
    pfData->getValue("BUS_TYPE", &ival);
    //printf ("bus_type = %d \n ", ival);

    if (ival != 4) {
      itmp ++;
      if ( temp*sin(rval) < 0) {
        printf("Powerflow bus%d mag = %f %fi\n", itmp, temp*cos(rval), temp*sin(rval));
      } else {
        printf("Powerflow bus%d mag = %f +%fi\n", itmp, temp*cos(rval), temp*sin(rval));
      }
    }
    int ngen = 0;
    pfData->getValue(GENERATOR_NUMBER, &ngen);
//    printf("number of gens = %d \n", ngen);
    if (pfData->getValue(GENERATOR_NUMBER, &ngen)) {
      int j;
      for (j=0; j<ngen; j++) {
        pfData->getValue("GENERATOR_PF_PGEN",&rval,j);
        dsData->setValue(GENERATOR_PG,rval,j);
        if (ngen >1) printf("save PGEN: %f\n", rval);
        pfData->getValue("GENERATOR_PF_QGEN",&rval,j);
        dsData->setValue(GENERATOR_QG,rval,j);
        if (ngen > 1) printf("save QGEN: %f\n", rval);
      }
    }
  }
}

// Calling program for the dynamis simulation applications

main(int argc, char **argv)
{
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);

  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

#ifdef USE_GOSS
  activemq::library::ActiveMQCPP::initializeLibrary();
#endif

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
    //multiply(*ybus, *volt_full, *INorton_full_chk);
    pf_app.write();
    pf_app.saveData();
   
/*    // setup and run dynamic simulation calculation
    boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork> ds_network(new
        gridpack::dynamic_simulation::DSFullNetwork(world));
    gridpack::dynamic_simulation::DSFullApp ds_app;
    int t_rn = timer->createCategory(" *** Dynamic Simulation: readNetwork ***");
    timer->start(t_rn);
    ds_app.readNetwork(ds_network,config);
    timer->stop(t_rn);
    int t_rg = timer->createCategory(" *** Dynamic Simulation: readGenerators ***");
    timer->start(t_rg);
    ds_app.readGenerators();
    timer->stop(t_rg);
    gridpack::utility::Configuration::CursorPtr cursor; cursor =
      config->getCursor("Configuration.Dynamic_simulation");
    std::vector<gridpack::dynamic_simulation::DSFullBranch::Event> faults;
    faults = ds_app.getFaults(cursor);
    int t_init = timer->createCategory(" *** Dynamic Simulation: initialization ***");
    timer->start(t_init);
    ds_app.initialize();
    timer->stop(t_init);
    int t_integration = timer->createCategory(" *** Dynamic Simulation: integration ***");
    timer->start(t_integration);
    ds_app.solve(faults[0]);
    timer->stop(t_integration);
    timer->stop(t_total);
    //timer->dump();*/
    
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
    std::vector<gridpack::dynamic_simulation::DSFullBranch::Event> faults;
    faults = ds_app.getFaults(cursor);

    // run dynamic simulation
    ds_app.setNetwork(ds_network, config);
    //ds_app.readNetwork(ds_network,config);
    printf("ds_app.readGenerators:\n");
    printf("initial number of buses: %d number of branches: %d\n",
        ds_network->numBuses(),ds_network->numBranches());
    ds_app.readGenerators();
    printf("number of buses: %d number of branches: %d\n",
        ds_network->numBuses(),ds_network->numBranches());
    printf("ds_app.initialize:\n");
    ds_app.initialize();
    printf("ds_app.setGenWatch:\n");
    ds_app.setGeneratorWatch();
    ds_app.setLoadWatch();
    ds_app.open("init_debug.values");
    ds_app.write("init_debug");
    ds_app.close();
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

