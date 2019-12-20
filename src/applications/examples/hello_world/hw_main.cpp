/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hw_main.cpp
 * @author Bruce Palmer
 * @date   2014-12-09 14:40:37 d3g096
 * 
 * @brief  
 */
// -------------------------------------------------------------

#include "gridpack/include/gridpack.hpp"
#include "hw_app.hpp"

// Calling program for the hello_world applications

int
main(int argc, char **argv)
{
  gridpack::Environment env(argc, argv);

  gridpack::hello_world::HWApp app;
  app.execute(argc, argv);

  return 0;
}
