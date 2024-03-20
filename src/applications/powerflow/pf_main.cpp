
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_main.cpp
 * @author Bruce Palmer
 * @date   2016-07-14 14:23:07 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"

const char* help = "GridPACK power flow application";

int main(int argc, char **argv)
{
  // Initialize libraries (parallel and math)
  gridpack::Environment env(argc,argv,help);

  if (1) {
    gridpack::utility::CoarseTimer *timer =
      gridpack::utility::CoarseTimer::instance();
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

    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config->getCursor("Configuration.Powerflow");
    bool useNonLinear = false;
    useNonLinear = cursor->get("UseNonLinear", useNonLinear);
    bool exportPSSE23 = false;
    std::string filename23;
    exportPSSE23 = cursor->get("exportPSSE_v23",&filename23);
    bool exportPSSE33 = false;
    std::string filename33;
    exportPSSE33 = cursor->get("exportPSSE_v33",&filename33);
    bool exportPSSE34 = false;
    std::string filename34;
    exportPSSE34 = cursor->get("exportPSSE_v34",&filename34);
    bool noPrint = false;
    cursor->get("suppressOutput",&noPrint);

    // setup and run powerflow calculation
    boost::shared_ptr<gridpack::powerflow::PFNetwork>
      pf_network(new gridpack::powerflow::PFNetwork(world));

    gridpack::powerflow::PFAppModule pf_app;
    if (noPrint) {
      pf_app.suppressOutput(noPrint);
    }
    pf_app.readNetwork(pf_network,config);
    pf_app.initialize();
    int ngen = pf_network->numGenerators();
    int nload = pf_network->numLoads();
    int nline = pf_network->numLines();
    if (world.rank() == 0) {
      std::cout<<"Number of generators in network: "<<ngen<<std::endl;
      std::cout<<"Number of loads in network:      "<<nload<<std::endl;
      std::cout<<"Number of lines in network:      "<<nline<<std::endl;
    }
    if (useNonLinear) {
      pf_app.nl_solve();
    } else {
      pf_app.solve();
      //pf_app.write();
    }
    pf_app.write();
    pf_app.saveData();
    if (exportPSSE23) {
      pf_app.exportPSSE23(filename23);
    }
    if (exportPSSE33) {
      pf_app.exportPSSE33(filename33);
    }
    if (exportPSSE34) {
      pf_app.exportPSSE34(filename34);
    }
    if (!noPrint) {
      timer ->dump();
    }
  }

  return 0;
}

