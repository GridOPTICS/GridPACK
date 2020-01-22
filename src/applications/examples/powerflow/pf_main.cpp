/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_main.cpp
 * @author Bruce Palmer
 * @date   2014-12-09 14:39:50 d3g096
 * 
 * @brief  
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "pf_app.hpp"

// Calling program for the powerflow application

int
main(int argc, char **argv)
{
  // Initialize parallel environment. This will start up
  // the underlying communication libraries
  gridpack::Environment env(argc,argv);

  // Initialize Math libraries
  gridpack::math::Initialize(&argc,&argv);

  // Create the power flow application and execute it
  gridpack::powerflow::PFApp app;
  app.execute(argc, argv);

  // Terminate Math libraries
  gridpack::math::Finalize();

  return 0;
}
