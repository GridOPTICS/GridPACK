// -------------------------------------------------------------
/**
 * @file   shuffle.cpp
 * @author William A. Perkins
 * @date Fri Mar 30 14:01:24 2012
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
// Last Change: Fri Mar 30 14:01:24 2012 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <ctime>
#include <iostream>
#include <vector>
#include <string>
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

#define VEC_LENGTH 5

struct Thing {
  int value;
  int orig_owner;
  int owner;
  int localidx;
  int claimed;
  std::vector<std::string> twister;

  /// Serialization method
  template <class Archive>
  void serialize(Archive &ar, const unsigned int)
  {
    ar & value;
    ar & orig_owner;
    ar & owner;
    ar & localidx;
    ar & claimed;
    ar & twister;
  }
};

std::ostream& 
operator<< (std::ostream& os, const Thing& t)
{
  return os << t.value;
}

typedef std::vector<Thing> ThingVector;

typedef std::vector<int> ivector;

  // mpi::all_gather(world, &mythings[0], mythings.size(), allthings);
  // if (world.rank() == 0) {
  //   std::cout << " src: ";
  //   for (ThingVector::iterator t = allthings.begin(); 
  //        t != allthings.end(); ++t) {
  //     std::cout << t->owner;
  //   }
  //   std::cout << std::endl;
  //   std::cout << "dest: ";
  //   for (ThingVector::iterator t = allthings.begin(); 
  //        t != allthings.end(); ++t) {
  //     std::cout << t->claimed;
  //   }
  //   std::cout << std::endl;
  // }
  

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
  std::vector<std::vector<std::string> > lstring;

  boost::mt19937 rengine(time(0) + world.rank());
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > 
    my_generator(rengine, boost::uniform_int<>(0, world.size() - 1));

  ThingVector mythings(num_per_proc);
  int gidx(num_per_proc*world.rank());
  int localidx(0);
  for (ThingVector::iterator t = mythings.begin(); 
       t != mythings.end(); ++t, ++gidx, ++localidx) {
    t->value = gidx;
    t->localidx = localidx;
    t->orig_owner = world.rank();
    t->owner = world.rank();
    t->claimed = my_generator();
    src.push_back(t->owner);
    dest.push_back(t->claimed);
    lidx.push_back(localidx);
    for (int i=0; i<VEC_LENGTH; i++) {
      char buf[128];
      sprintf(buf,"gidx: %d localidx: %d orig_owner: %d i: %d",
          gidx,localidx,world.rank(),i);
      t->twister.push_back(buf);
    }
  }

  std::cout << world.rank() << ": ";
  std::copy(mythings.begin(), mythings.end(),
            std::ostream_iterator<Thing>(std::cout, ", "));
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
  std::list<Thing> results;
  std::vector<mpi::request> reqs;

  for (; s != allsrc.end(); ++s, ++d, ++msgid, ++i) {
    if (world.rank() == *s) {
      localidx = *i;
      // std::cout << world.rank() << ": " << msgid << ": "
      //           << *s << " --> " << *d 
      //           << std::endl;
      reqs.push_back(world.isend(*d, msgid, &mythings[localidx], 1));
    }
    if (world.rank() == *d) {
      Thing t;
      results.push_back(t);
      reqs.push_back(world.irecv(*s, msgid, &(results.back()), 1));
      // std::cout << world.rank() << ": " << msgid << ": "
      //           << *d << " <-- " << *s
      //           << std::endl;
    } 
  }
  mpi::wait_all(reqs.begin(), reqs.end());
  world.barrier();

  std::cout << world.rank() << ": ";
  std::copy(results.begin(), results.end(),
           std::ostream_iterator<Thing>(std::cout, ", "));
  std::cout << std::endl;

  world.barrier();
  return 0;
  
}
