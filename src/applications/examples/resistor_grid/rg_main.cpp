/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   rg_main.cpp
 * @author Bruce Palmer
 * @date   2014-12-09 14:41:19 d3g096
 * 
 * @brief  
 */
// -------------------------------------------------------------

#include "gridpack/include/gridpack.hpp"
#include "rg_app.hpp"

// Calling program for the resistor_grid applications

int
main(int argc, char **argv)
{
  gridpack::Environment env(argc, argv);

  gridpack::resistor_grid::RGApp app;
  app.execute(argc, argv);

  return 0;
}
