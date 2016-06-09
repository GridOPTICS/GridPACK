// -------------------------------------------------------------
// file: shuffler_test_1.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created June 24, 2013 by William A. Perkins
// Last Change: 2013-12-16 12:26:37 d3g096
// -------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

boost::random::mt19937 gen;

#include "shuffler.hpp"
#include "printit.hpp"

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  const int local_size(8);
  size_t global_size(0);

  gen.seed(time(NULL)+world.rank());

  std::vector<int> src(local_size, world.rank());
  std::vector<std::string> mine;
  std::vector<int> dest;
  mine.resize(local_size);
  dest.resize(local_size);
  for (int i = 0; i < mine.size(); ++i) {
    char c('A');
    c += (i) % 26;
    size_t l(i % local_size + 1);
    mine[i] = std::string(5*l, c);
    dest[i] = i % world.size();
  }

  Shuffler<std::string, int, false> shuffle;

  printit(world, mine, "Before: ");
  shuffle(world, mine, dest);
  printit(world, mine, "After: ");
  
  dest.clear();
  dest.resize(mine.size(), 0);

  shuffle(world, mine, dest);
  printit(world, mine, "Back: ");

  dest.clear();
  dest.resize(mine.size(), 0);
  for (int i = 0; i < mine.size(); ++i) {
    dest[i] = i % world.size();
  }

  shuffle(world, mine, dest);
  printit(world, mine, "Again: ");

  std::cout << world.rank() << ": exiting " << std::endl;
  world.barrier();
  return 0;
}

