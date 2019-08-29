/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   KalmanDS_main.cpp
 * @author Da Meng and Yousu Chen 
 * @date   1/06/2015
 *
 * @brief
 *
 * Modified by Xinya Li, July 2015
 */ 
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "gridpack/applications/modules/kalman_ds/kds_app_module.hpp"

// Calling program for the state estimation applications

int main(int argc, char **argv)
{
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);
  // Intialize Math libraries
  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  gridpack::math::Initialize(&argc,&argv);

  {
    gridpack::utility::CoarseTimer *timer =
      gridpack::utility::CoarseTimer::instance();
    int t_Total = timer->createCategory("App:Total");
    int t_PF = timer->createCategory("PF: Total");
    int t_KF = timer->createCategory("KF: Time Loop");
    int t_In = timer->createCategory("KF: Input and Initialization");
    timer->start(t_Total);
    timer->start(t_PF);

    gridpack::parallel::Communicator comm;
    // Initialize Kalman filter calculation by first running a powerflow
    // simulation
    boost::shared_ptr<gridpack::powerflow::PFNetwork>
      pf_network(new gridpack::powerflow::PFNetwork(comm));

    // Read configuration file
    gridpack::utility::Configuration *config =
      gridpack::utility::Configuration::configuration();
    if (argc >= 2 && argv[1] != NULL) {
      char inputfile[256];
      sprintf(inputfile,"%s",argv[1]);
      config->open(inputfile,comm);
    } else {
      config->open("input.xml",comm);
    }

    // run powerflow calculation and save data to data collection objects
    gridpack::powerflow::PFAppModule pf_app;
    pf_app.readNetwork(pf_network,config);
    pf_app.initialize();
    pf_app.solve();
    pf_app.write();
    pf_app.saveData();
    timer->stop(t_PF);

    boost::shared_ptr<gridpack::kalman_filter::KalmanNetwork>
      kds_network(new gridpack::kalman_filter::KalmanNetwork(comm));
    pf_network->clone<gridpack::kalman_filter::KalmanBus,
          gridpack::kalman_filter::KalmanBranch>(kds_network);
    gridpack::kalman_filter::KalmanApp kds_app;
    kds_app.setNetwork(kds_network, config);
    kds_app.initialize();
    kds_app.solve();
    timer->stop(t_Total);
    timer->dump();
  }

  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
}

