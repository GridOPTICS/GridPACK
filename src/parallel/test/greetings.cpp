// -------------------------------------------------------------
/**
 * @file   greetings.cpp
 * @author William A. Perkins
 * @date   2013-05-07 09:58:11 d3g096
 * 
 * @brief  A simple test of the GridPACK parallel environment
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created May  6, 2013 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------

#include <iostream>
#include "gridpack/parallel/parallel.hpp"

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  gridpack::parallel::Communicator world;

  std::cout << "I am process " << world.rank() 
            << " of " << world.size()
            << "." << std::endl;
  
  return 0;
}

