/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_main.cpp
 * @author Bruce Palmer
 * @date   2016-07-14 14:23:28 d3g096
 *
 * @brief
 * @update Yousu Chen
 *         Adding functions of applying virtual measurements, enhancing voltage 
 *         constraings, and performing bad data detection
 * @date   2025-03-05
 *         Adding try/catch function to enable more robust state estimation
 * @date   2025-04-20
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "timer/coarse_timer.hpp"

// Calling program for the state estimation application

int
main(int argc, char **argv)
{
  int return_code = 0;
  
  try {
    // Initialize MPI libraries
    int ierr = MPI_Init(&argc, &argv);

    GA_Initialize();
    int stack = 200000, heap = 200000;
    MA_init(C_DBL, stack, heap);

    // Intialize Math libraries
    gridpack::math::Initialize(&argc,&argv);

    try {
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

      // setup and run state estimation calculation
      boost::shared_ptr<gridpack::state_estimation::SENetwork>
        se_network(new gridpack::state_estimation::SENetwork(world));

      gridpack::state_estimation::SEAppModule se_app;
      
      try {
        se_app.readNetwork(se_network,config);
      } catch (const std::exception& e) {
        if (world.rank() == 0) {
          printf("ERROR: Exception during network reading: %s\n", e.what());
        }
        return_code = 1;
        throw; // Re-throw to outer handler
      } catch (...) {
        if (world.rank() == 0) {
          printf("ERROR: Unknown exception during network reading\n");
        }
        return_code = 1;
        throw; // Re-throw to outer handler
      }
      
      try {
        se_app.initialize();
      } catch (const std::exception& e) {
        if (world.rank() == 0) {
          printf("ERROR: Exception during initialization: %s\n", e.what());
        }
        return_code = 1;
        throw; // Re-throw to outer handler
      } catch (...) {
        if (world.rank() == 0) {
          printf("ERROR: Unknown exception during initialization\n");
        }
        return_code = 1;
        throw; // Re-throw to outer handler
      }
      
      try {
        se_app.readMeasurements();
      } catch (const std::exception& e) {
        if (world.rank() == 0) {
          printf("ERROR: Exception during measurement reading: %s\n", e.what());
        }
        return_code = 1;
        throw; // Re-throw to outer handler
      } catch (...) {
        if (world.rank() == 0) {
          printf("ERROR: Unknown exception during measurement reading\n");
        }
        return_code = 1;
        throw; // Re-throw to outer handler
      }
      
      try {
        printf("DEBUG: About to call se_app.solve()\n");
        fflush(stdout);
        se_app.solve();
        printf("DEBUG: se_app.solve() completed, about to call saveData()\n");
        fflush(stdout);
        // Save data to data collection objects
        se_app.saveData();
        printf("DEBUG: se_app.saveData() completed\n");
        fflush(stdout);
      } catch (const std::exception& e) {
        if (world.rank() == 0) {
          printf("ERROR: Exception caught during state estimation: %s\n", e.what());
          printf("Continuing to write available results...\n");
        }
        return_code = 1;
      } catch (...) {
        if (world.rank() == 0) {
          printf("ERROR: Unknown exception caught during state estimation\n");
          printf("Continuing to write available results...\n");
        }
        return_code = 1;
      }
      
      printf("DEBUG: About to call se_app.write()\n");
      fflush(stdout);
      try {
        // Write results regardless of convergence status
        se_app.write();
        printf("DEBUG: se_app.write() completed\n");
        fflush(stdout);
      } catch (const std::exception& e) {
        if (world.rank() == 0) {
          printf("ERROR: Exception during result writing: %s\n", e.what());
        }
        return_code = 1;
      } catch (...) {
        if (world.rank() == 0) {
          printf("ERROR: Unknown exception during result writing\n");
        }
        return_code = 1;
      }
      
      printf("DEBUG: About to dump timing information\n");
      fflush(stdout);
      // Dump timing information
      if (world.rank() == 0) {
        gridpack::utility::CoarseTimer::instance()->dump();
      }
      printf("DEBUG: Timing dump completed\n");
      fflush(stdout);
      
      // Report convergence status
      if (!se_app.hasConverged()) {
        if (world.rank() == 0) {
          printf("WARNING: State estimation did not fully converge, but results were written anyway.\n");
        }
      }
      printf("DEBUG: Convergence check completed\n");
      fflush(stdout);
    } catch (const std::exception& e) {
      if (argc >= 1 && argv[0] != NULL) {
        printf("%s: ERROR: Unhandled exception: %s\n", argv[0], e.what());
      } else {
        printf("ERROR: Unhandled exception: %s\n", e.what());
      }
      return_code = 1;
    } catch (...) {
      if (argc >= 1 && argv[0] != NULL) {
        printf("%s: ERROR: Unknown unhandled exception\n", argv[0]);
      } else {
        printf("ERROR: Unknown unhandled exception\n");
      }
      return_code = 1;
    }

    // Clean up regardless of exceptions
    GA_Terminate();
    gridpack::math::Finalize();
    MPI_Finalize();
    
  } catch (...) {
    // Last-resort catch for any exceptions during initialization/finalization
    if (argc >= 1 && argv[0] != NULL) {
      printf("%s: FATAL ERROR: Exception during program initialization or finalization\n", argv[0]);
    } else {
      printf("FATAL ERROR: Exception during program initialization or finalization\n");
    }
    return_code = 2;
  }
  
  return return_code;
}

