// -------------------------------------------------------------
// file: shuffler_test_1.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created June 24, 2013 by William A. Perkins
// Last Change: 2013-12-16 07:57:00 d3g096
// -------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

boost::random::mt19937 gen;

#include "shuffler.hpp"
#include "printit.hpp"

// -------------------------------------------------------------
// redistribute_randomly
// -------------------------------------------------------------
template <typename Thing>
struct RandomShuffler {
  typedef std::vector<Thing> Vector;
  // typedef std::vector<Thing>::difference_type Index;
  typedef int Index;
  typedef std::vector<Index> IndexVector;

  IndexVector dest;
  
  void
  operator() (const boost::mpi::communicator& comm,
              boost::random::mt19937& gen,
              Vector& mine)
  {
    dest.clear();
    if (!mine.empty()) {
      dest.resize(mine.size());
      boost::random::uniform_int_distribution<> dist(0, comm.size()-1);
      for (int i = 0; i < mine.size(); ++i) {
        dest[i] = dist(gen);
      }
    }
    Shuffler<Thing, int, false> shuffle;
    shuffle(comm, mine, dest);
  }
};

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
  std::vector<int> mine(local_size);
  for (int i = 0; i < mine.size(); ++i) {
    mine[i] = world.rank()*local_size + i;
  }

  RandomShuffler<int> shuffle;

  printit(world, mine, "Before: ");
  shuffle(world, gen, mine);
  printit(world, mine, "After: ");

  std::vector<int> dest(shuffle.dest);
  Shuffler<int> reshuffle;
  reshuffle(world, src, dest);
  reshuffle(world, mine, src);

  printit(world, mine, "Return: ");

  return 0;
}

