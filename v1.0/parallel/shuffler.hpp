/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   shuffler.hpp
 * @author William A. Perkins
 * @date   2013-12-16 14:23:26 d3g096
 * 
 * @brief  A thing to redistribute a vector of things over several processors 
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
#include <utility>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>


#ifndef _shuffler_hpp_
#define _shuffler_hpp_

// -------------------------------------------------------------
//  class Shuffler
// -------------------------------------------------------------
/// A functor to redistribute things in vectors on multiple processes
/**
 *
 * Each process starts with a (possibly empty) vector of things and a
 * vector of equal size that containing a destination process for each
 * thing.  After execution, each process will contain a vector of the
 * things assigned to it.
 *
 * This uses blocking send/receive.  
 *
 * The things redistributed must be copy constructable and serializable.  
 * 
 */

template <typename Thing, typename I = int>
struct Shuffler {

  // Thing must be copyable, at a relatively low cost, and
  // serializable -- should test for that

  typedef std::vector<Thing> ThingVector;
  typedef I Index;
  typedef std::vector<Index> IndexVector;
  
  /// Redistribute and get the Things assigned to the local process
  void operator()(const boost::mpi::communicator& comm, 
                  ThingVector& locthings, const IndexVector& destproc)
  {
    BOOST_ASSERT(locthings.size() == destproc.size());

    size_t nthings(0);
    all_reduce(comm, locthings.size(), nthings, std::plus<size_t>());
    if (nthings <= 0) return;

    // save the original list of local things 

    ThingVector tvect; 
    tvect.reserve(locthings.size());
    std::copy(locthings.begin(), locthings.end(),
              std::back_inserter(tvect));
    locthings.clear();

    // Work on each processors list of things separately

    std::vector<ThingVector> tosend(comm.size());

    // all processes go through the destinations and makes a vector to
    // send to each of the other processes

    size_t locidx(0);

    for (typename IndexVector::const_iterator dest = destproc.begin(); 
         dest != destproc.end(); ++dest) {
      if (*dest == comm.rank()) {
        locthings.push_back(tvect[locidx]);
      } else {
        tosend[*dest].push_back(tvect[locidx]);
      }
      locidx += 1;
    }

    for (int src = 0; src < comm.size(); ++src) {
      ThingVector tmp;
      if (comm.rank() == src) {
        scatter(comm, tosend, tmp, src);
      } else {
        scatter(comm, tmp, src);
      }
      std::copy(tmp.begin(), tmp.end(), std::back_inserter(locthings));
    }
  }

};



#endif
