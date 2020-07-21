
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
    bool exportPSSE = false;
    std::string filename;
    exportPSSE = cursor->get("exportPSSE_v33",&filename);
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
    if (useNonLinear) {
      pf_app.nl_solve();
    } else {
      pf_app.solve();
      //pf_app.write();
    }
    pf_app.write();
    pf_app.saveData();
    if (exportPSSE) {
      pf_app.exportPSSE33(filename);
    }
    if (!noPrint) {
      timer ->dump();
    }
  }

  return 0;
}

