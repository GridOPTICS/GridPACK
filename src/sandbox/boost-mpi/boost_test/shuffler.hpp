// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   shuffler.hpp
 * @author William A. Perkins
 * @date   2015-12-18 08:01:12 d3g096
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

#ifndef _shuffler_hpp_
#define _shuffler_hpp_

#include <iostream>
#include <vector>
#include <list>
#include <utility>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>


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

template <typename Thing, typename I = int, bool doBlocking = true>
struct Shuffler {

  // Thing must be copyable and serializable and have a constructor
  // without arguments -- should test for that

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

    // save the original list of local things 

    ThingVector tmpthings(locthings); 
    locthings.clear();

    // Work on each processors list of things separately

    int msgid(0);               // unique MPI message id
    iPairVector srcdest;
    std::list<boost::mpi::request> reqs;
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
           pd != srcdest.end(); ++pd, ++msgid) {
        int src(pd->first);
        int dest(pd->second);
        if (comm.rank() == src) {
          if (src == dest) {
            locthings.push_back(tmpthings[locidx]);
            // std::cout << src << ": kept " << locthings.back() << " @ " 
            //           << locthings.size() - 1 << std::endl;
          } else {
            if (doBlocking) {
              comm.send(dest, msgid, tmpthings[locidx]);
            } else {
              reqs.push_back(comm.isend(dest, msgid, &tmpthings[locidx], 1));
            }
            // std::cout << src << ": msg " << msgid << ": sent " << 
            //   tmpthings[locidx] << std::endl;
          }
          locidx += 1;
        } else if (comm.rank() == dest) {
          Thing bogus;
          locthings.push_back(bogus);
          if (doBlocking) {
            comm.recv(src, msgid, locthings.back());
          } else {
            reqs.push_back(comm.irecv(src, msgid, &(locthings.back()), 1));
          }
          // std::cout << dest << ": msg " << msgid 
          //           << ": got @ " << locthings.size() - 1 << " \"" 
          //           << locthings.back() << "\"" << std::endl;
        }
        comm.barrier();
      }

      if (!doBlocking) {
        // std::cout << comm.rank() << ": "
        //           << reqs.size() << " non-blocking messages created"
        //           << " (src = " << p << ")"
        //           << std::endl;
        boost::mpi::wait_all(reqs.begin(), reqs.end());
        reqs.clear();
      }
      
      comm.barrier();
    }

    tmpthings.clear();
  }

};



#endif
