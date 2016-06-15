// -------------------------------------------------------------
/**
 * @file   shuffle.cpp
 * @author William A. Perkins
 * @date Fri Mar 30 14:01:03 2012
 * 
 * @brief  A program to test point-to-point communication
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 29, 2012 by William A. Perkins
// Last Change: Fri Mar 30 14:01:03 2012 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <ctime>
#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <cstdlib>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/nonblocking.hpp>

namespace mpi = boost::mpi;

#include <boost/serialization/vector.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>

typedef std::vector<int> ivector;

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  static const int num_per_proc(4);

  ivector src, allsrc;
  ivector dest, alldest;
  ivector lidx, alllidx;

  boost::mt19937 rengine(time(0)+world.rank());
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > 
    my_generator(rengine, boost::uniform_int<>(0, world.size() - 1));

  ivector myvalues(num_per_proc, -1);
  int gidx(num_per_proc*world.rank());
  int localidx(0);

  for (ivector::iterator t = myvalues.begin(); 
       t != myvalues.end(); ++t, ++gidx, ++localidx) {
    *t = gidx;
    src.push_back(world.rank());
    dest.push_back(my_generator());
    lidx.push_back(localidx);
  }

  std::cout << world.rank() << ": ";
  std::copy(myvalues.begin(), myvalues.end(),
            std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;

  mpi::all_gather(world, &src[0], src.size(), allsrc);
  mpi::all_gather(world, &dest[0], dest.size(), alldest);
  mpi::all_gather(world, &lidx[0], lidx.size(), alllidx);

  if (world.rank() == 0) {
    std::cout << " src: ";
    std::copy(allsrc.begin(), allsrc.end(), 
              std::ostream_iterator<int>(std::cout));
    std::cout << std::endl;
    std::cout << "dest: ";
    std::copy(alldest.begin(), alldest.end(), 
              std::ostream_iterator<int>(std::cout));
    std::cout << std::endl;
    std::cout << "lidx: ";
    std::copy(alllidx.begin(), alllidx.end(), 
              std::ostream_iterator<int>(std::cout));
    std::cout << std::endl;
  }

  world.barrier();
  int msgid(0);
  ivector::iterator s(allsrc.begin());
  ivector::iterator d(alldest.begin());
  ivector::iterator i(alllidx.begin());
  std::list<int> results;
  std::vector<mpi::request> reqs;

  for (; s != allsrc.end(); ++s, ++d, ++msgid, ++i) {
    if (world.rank() == *s) {
      localidx = *i;
      std::cout << world.rank() << ": " << msgid << ": "
                << *s << " --> " << *d 
                << std::endl;
      reqs.push_back(world.isend(*d, msgid, myvalues[localidx]));
    }
    if (world.rank() == *d) {
      results.push_back(-1);
      reqs.push_back(world.irecv(*s, msgid, results.back()));
      std::cout << world.rank() << ": " << msgid << ": "
                << *d << " <-- " << *s
                << std::endl;
    } 
  }
  mpi::wait_all(reqs.begin(), reqs.end());
  world.barrier();

  std::cout << world.rank() << ": ";
  std::copy(results.begin(), results.end(),
            std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;

  world.barrier();
  return 0;
  
}
