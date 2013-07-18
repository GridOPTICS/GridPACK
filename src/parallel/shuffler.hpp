// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   shuffler.hpp
 * @author William A. Perkins
 * @date   2013-07-15 10:28:44 d3g096
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
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/nonblocking.hpp>
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
 * The things redistributed must be copy constructable and serializable.  
 * 
 */

template <typename Thing, typename I = int>
struct Shuffler {

  // Thing must be copyable and serializable -- should test for that

  typedef std::vector<Thing> ThingVector;
  typedef I Index;
  typedef std::vector<Index> IndexVector;
  typedef std::pair<Index, Index> iPair;
  typedef std::vector<iPair> iPairVector;
  
  /// Redistribute and get the Things assigned to the local process
  void operator()(const boost::mpi::communicator& comm, 
                  ThingVector& locthings, const IndexVector& destproc)
  {
    BOOST_ASSERT(locthings.size() == destproc.size());

    size_t nthings(0);
    all_reduce(comm, locthings.size(), nthings, std::plus<size_t>());
    if (nthings <= 0) return;

    ThingVector tvect(locthings); 
    
    locthings.clear();

    int msgid(0);
    std::list<boost::mpi::request> reqs;
    
    // Work on each processors list of things separately
    iPairVector srcdest;
    for (int p = 0; p < comm.size(); ++p) {
      
      // build a set of source-destination pairs for each thing in p's
      // local list

      srcdest.clear();
      if (comm.rank() == p) {
        for (typename IndexVector::const_iterator i = destproc.begin(); 
             i != destproc.end(); ++i) {
          iPair apair(p, *i);
          srcdest.push_back(apair);
        }
      }

      // all processes get the source-destination pairs

      broadcast(comm, srcdest, p);

      // all processes go through the source-destination pairs and
      // build a message request list 

      size_t locidx(0);
      for (typename iPairVector::iterator pd = srcdest.begin(); 
           pd != srcdest.end(); ++pd, ++msgid, ++locidx) {
        int src(pd->first);
        int dest(pd->second);
        if (src == dest && comm.rank() == src) {
            locthings.push_back(tvect[locidx]);
        } else if (comm.rank() == src) {
          reqs.push_back(comm.isend(dest, msgid, &tvect[locidx], 1));
        } else if (comm.rank() == dest) {
          // push an empty thing on the local list so there is some
          // place to recieve it
          Thing bogus;
          locthings.push_back(bogus);
          reqs.push_back(comm.irecv(src, msgid, &(locthings.back()), 1));
        }
        comm.barrier();
      }

      // execute the message requests 
      
      boost::mpi::wait_all(reqs.begin(), reqs.end());

      // this needs to be here in case a processor was not involved in
      // any transfer

      comm.barrier();
    }
  }

};



#endif
