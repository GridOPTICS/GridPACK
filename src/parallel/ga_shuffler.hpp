// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   ga_shuffler.hpp
 * @author William A. Perkins
 * @date   2014-03-05 10:41:05 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _ga_shuffler_hpp_
#define _ga_shuffler_hpp_

#include <string>
#include <cstring>
#include <sstream>
#include <vector>

#include <boost/assert.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <ga++.h>

#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/parallel/distributed.hpp>

namespace gridpack {
namespace parallel {

// -------------------------------------------------------------
//  class gaShuffler
// -------------------------------------------------------------
template <typename Thing, typename Index = int> 
class gaShuffler 
  : public Distributed,
    private utility::Uncopyable
{
public:

  // Thing must be copyable, at a relatively low cost, and
  // serializable -- should test for that

  typedef std::vector<Thing> ThingVector;
  typedef std::vector<Index> IndexVector;
  
  /// Default constructor.
  gaShuffler(const Communicator& comm)
    : Distributed(comm), utility::Uncopyable(),
      p_gaBuffers(processor_size()),
      p_gaBufferSize(processor_size())
  {}

  /// Destructor
  ~gaShuffler(void)
  {}

  /// Move some things to different processors
  void operator() (ThingVector& locthings, const IndexVector& destproc)
  {
    BOOST_ASSERT(locthings.size() == destproc.size());

    int me(processor_rank());
    int nproc(processor_size());

    if (nproc <= 1) return;

    size_t nthings(0);
    boost::mpi::all_reduce(this->communicator(), locthings.size(), nthings, 
                           std::plus<size_t>());
    if (nthings <= 0) return;

    int theGAgroup(communicator().getGroup());
    int oldGAgroup = GA_Pgroup_get_default();
    GA_Pgroup_set_default(theGAgroup);

    // save the original list of local things 

    ThingVector tvect(locthings); 
    locthings.clear();

    // Work on each processors list of things separately

    std::vector< ThingVector > tosend(processor_size());

    // all processes go through the destinations and makes a vector to
    // send to each of the other processes

    typename ThingVector::iterator t(tvect.begin());
    typename IndexVector::const_iterator dest(destproc.begin());
    for (; dest != destproc.end(); ++dest, ++t) {
      if (*dest == processor_rank()) {
        locthings.push_back(*t);
      } else {
        tosend[*dest].push_back(*t);
      }
    }

    for (int p = 0; p < nproc; ++p) {
      if (p != me) {
        std::ostringstream ostr(std::ios::binary);
        boost::archive::binary_oarchive oarch(ostr);
        oarch & tosend[p];
        int buflen = MA_sizeof(MT_CHAR, ostr.str().length(), p_gaType);
        p_gaBuffers[p].resize(buflen);
        memcpy(&(p_gaBuffers[p][0]), 
               ostr.str().c_str(), 
               ostr.str().length()*sizeof(char));
        p_gaBufferSize[p] = p_gaBuffers[p].size();
      }
    }

    size_t lmaxsize(0), thesize;
    lmaxsize = *std::max_element(p_gaBufferSize.begin(),
                                 p_gaBufferSize.end());
    boost::mpi::all_reduce(communicator(), lmaxsize, thesize, 
                           boost::mpi::maximum<size_t>());

    int dims[2] = { nproc, thesize + 1 };
    int lo[2], hi[2], ld[2] = {1, 1};
    boost::scoped_ptr<GA::GlobalArray> 
      ga(new GA::GlobalArray(p_gaType, 2, dims, "serialized data", NULL)),
      gasize(new GA::GlobalArray(MT_C_INT, 1, dims, "serialized data sizes", NULL));

    for (int src = 0; src < nproc; ++src) {
      if (me == src) {
        lo[0] = 0, hi[0] = nproc - 1;
        gasize->put(&lo[0], &hi[0], &(p_gaBufferSize[0]), &ld[0]);
        for (int p = 0; p < nproc; ++p) {
          if (p != me) {
            lo[0] = hi[0] = p;
            lo[1] = 0; hi[1] = p_gaBufferSize[p] - 1;
            ga->put(&lo[0], &hi[0], &(p_gaBuffers[p][0]), &ld[0]);
          }
        }
      }
      communicator().sync();
      if (me != src) {
        int mybufsize;
        lo[0] = hi[0] = me;
        gasize->get(&lo[0], &hi[0], &mybufsize, &ld[0]);
        lo[1] = 0; hi[1] = mybufsize - 1;
        p_gaBuffers[me].clear();
        p_gaBuffers[me].resize(mybufsize);
        ga->get(&lo[0], &hi[0], &(p_gaBuffers[me][0]), &ld[1]);
        int buflen(p_gaBuffers[me].size());
        std::string s((char *)(&(p_gaBuffers[me][0])), 
                      buflen*sizeof(gaBufferType::value_type));
        p_gaBuffers[me].clear();
        std::istringstream is(s, std::ios::binary);
        boost::archive::binary_iarchive iarch(is);
        ThingVector tmp;
        iarch & tmp;
        std::copy(tmp.begin(), tmp.end(), std::back_inserter(locthings));
      }
      communicator().sync();
    }

    ga.reset();
    gasize.reset();

    p_reset();

    GA_Pgroup_set_default(oldGAgroup);
  }

private:

  /// The GA type to use
  static const int p_gaType = MT_C_INT;

  /// Buffer type compatible with GA type
  typedef std::vector<int> gaBufferType;

  /// GA compatible buffers for each processor
  std::vector<gaBufferType> p_gaBuffers;

  /// The actual sizes of the entries in p_gaBuffers
  std::vector<int> p_gaBufferSize;

  /// Reset this instance
  void p_reset(void)
  {
    for (size_t p = 0; p < p_gaBuffers.size(); ++p) {
      p_gaBuffers[p].clear();
      p_gaBufferSize[p] = 0;
    }
  }
};

} // namespace gridpack
} // namespace parallel

#endif
