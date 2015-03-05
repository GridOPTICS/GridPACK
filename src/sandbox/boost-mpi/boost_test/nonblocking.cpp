// -------------------------------------------------------------
// file: nonblocking.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December  5, 2013 by William A. Perkins
// Last Change: 2013-12-05 12:30:40 d3g096
// -------------------------------------------------------------


// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
#include <boost/mpi.hpp>
#include <iostream>
#include <string>
#include <list>
#include <boost/serialization/string.hpp>

namespace mpi = boost::mpi;

int 
main(int argc, char* argv[])
{
  mpi::environment env(argc, argv);
  mpi::communicator world;
  
  std::list<mpi::request> reqs;
  std::string msg, out_msg = "Hello";
  if (world.rank() == 0) {
    out_msg = "Hello big";
    reqs.push_back(world.isend(1, 0, out_msg));
    reqs.push_back(world.irecv(1, 1, msg));
  } else {
    out_msg = "world";
    reqs.push_back(world.isend(0, 1, out_msg));
    reqs.push_back(world.irecv(0, 0, msg));
  }
  mpi::wait_all(reqs.begin(), reqs.end());
  std::cout << world.rank() << ": " << msg << std::endl;
  world.barrier();
  return 0;
}
