// -------------------------------------------------------------
/**
 * @file   greetings.cpp
 * @author William A. Perkins
 * @date Thu Mar 29 10:29:33 2012
 * 
 * @brief  A simple boost::mpi test program
 * 
 * 
 * This example come directly from the boost::mpi tutorial.  Each
 * process just prints a message.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 28, 2012 by William A. Perkins
// Last Change: Thu Mar 29 10:29:33 2012 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <iostream>

namespace mpi = boost::mpi;

int main(int argc, char* argv[])
{
  mpi::environment env(argc, argv);
  mpi::communicator world;
  std::cout << "I am process " << world.rank() << " of " << world.size()
            << "." << std::endl;
  return 0;
}
