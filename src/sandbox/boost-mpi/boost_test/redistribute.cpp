// -------------------------------------------------------------
/**
 * @file   redistribute.cpp
 * @author William A. Perkins
 * @date   2015-12-18 08:01:49 d3g096
 * 
 * @brief  Explore ways to redistribute vectors over processors
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created June 24, 2013 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include "shuffler.hpp"


// -------------------------------------------------------------
// redistribute_evenly
// -------------------------------------------------------------
template <typename T>
void
redistribute_evenly(const boost::mpi::communicator& comm,
                    std::vector<T>& locthings)
{
  std::vector<size_t> oldperproc(comm.size());
  std::vector<size_t> newperproc(comm.size(), 0);
  if (comm.rank() == 0) {
    gather(comm, locthings.size(), oldperproc, 0);
    size_t allsize = std::accumulate(oldperproc.begin(), oldperproc.end(), 0);
    for (size_t i = 0; i < allsize; ++i) {
      size_t p(i % comm.size());
      newperproc[p] += 1;
    }
    std::cout << "Old: ";
    std::copy(oldperproc.begin(), oldperproc.end(),
              std::ostream_iterator<size_t>(std::cout, ","));
    std::cout << std::endl;
    std::cout << "New: ";
    std::copy(newperproc.begin(), newperproc.end(),
              std::ostream_iterator<size_t>(std::cout, ","));
    std::cout << std::endl;
  } else {
    gather(comm, locthings.size(), 0);
  }
  broadcast(comm, newperproc, 0);

  
}

// -------------------------------------------------------------
// printit
// -------------------------------------------------------------
template <typename T> 
void
printit(const boost::mpi::communicator& comm, const std::vector<T> things)
{
  for (int p = 0; p < comm.size(); ++p) {
    if (comm.rank() == p) {
      std::cout << p << ": ";
      std::copy(things.begin(), things.end(),
                std::ostream_iterator<T>(std::cout, ","));
      std::cout << std::endl;
      std::cout.flush();
    }
    comm.barrier();
  }
}                                      

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  const size_t local_size(8);
  const size_t global_size(world.size()*local_size - 1);
  
  std::vector<int> mine;
  if (world.rank() == 0) {
    mine.reserve(global_size);
    std::copy(boost::counting_iterator<int>(0),
              boost::counting_iterator<int>(global_size),
              std::back_inserter(mine));
  }

  printit(world, mine);

  redistribute_evenly(world, mine);

  printit(world, mine);
  
  return 0;
}

