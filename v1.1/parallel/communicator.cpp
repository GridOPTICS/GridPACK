// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   communicator.cpp
 * @author William A. Perkins
 * @date   2014-02-13 09:22:41 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <ga++.h>
#include "gridpack/utilities/uncopyable.hpp"
#include "communicator.hpp"

namespace gridpack {
namespace parallel {

// -------------------------------------------------------------
//  class CommunicatorPrivate
// -------------------------------------------------------------
class CommunicatorPrivate {
public:

  /// Default constructor.
  CommunicatorPrivate(void)
    : p_handle(GA_Pgroup_get_world())
  { }

  /// Construct on the processes in the specified communicators
  CommunicatorPrivate(const boost::mpi::communicator& comm)
    : p_handle()
  {
    int me(comm.rank()), nprocs(comm.size());
    int gaWrld = GA_Pgroup_get_world();
    int defGrp = GA_Pgroup_get_default();
    GA_Pgroup_set_default(gaWrld);
    std::vector<int> gaSrc(nprocs, 0), gaDest(nprocs, 0);
    gaSrc[me] = GA_Nodeid();
    boost::mpi::all_reduce(comm, &gaSrc[0], nprocs, &gaDest[0], std::plus<int>());
    p_handle = GA_Pgroup_create(&gaDest[0],nprocs);
    GA_Pgroup_sync(p_handle);
    GA_Pgroup_set_default(defGrp);
  }

  /// Destructor
  ~CommunicatorPrivate(void)
  { 
    if (p_handle != GA_Pgroup_get_world()){
      GA_Pgroup_destroy(p_handle);
    }
  }

  /// Get GA process group handle
  const int& handle(void) const
  {
    return p_handle;
  }

private:

  /// The GA process group handle
  int p_handle;
  
};


// -------------------------------------------------------------
//  class Communicator
// -------------------------------------------------------------

// -------------------------------------------------------------
// Communicator:: constructors / destructor
// -------------------------------------------------------------
Communicator::Communicator(void)
  : p_comm(), p_private(new CommunicatorPrivate())
{
  
}

Communicator::Communicator(const Communicator& old)
  : p_comm(old.p_comm), p_private(old.p_private)
{
}

Communicator::Communicator(const boost::mpi::communicator& comm)
  : p_comm(comm), p_private(new CommunicatorPrivate(p_comm))
{
}


Communicator::Communicator(const MPI_Comm& comm)
  : p_comm(comm, boost::mpi::comm_duplicate),
    p_private(new CommunicatorPrivate(p_comm))
{
}


Communicator::~Communicator(void)
{
}

// -------------------------------------------------------------
// Communicator::swap
// -------------------------------------------------------------
void
Communicator::swap(Communicator& other)
{
  std::swap(p_comm, other.p_comm);
  std::swap(p_private, other.p_private);
}

// -------------------------------------------------------------
// Communicator::operator=
// -------------------------------------------------------------
Communicator & 
Communicator::operator= (const Communicator & rhs) 
{
  Communicator temp(rhs);
  this->swap(temp);
  return (*this);
}

// -------------------------------------------------------------
// Communicator::worldRank
// -------------------------------------------------------------
int 
Communicator::worldRank(void) const
{
  return GA::SERVICES.nodeid();
}

// -------------------------------------------------------------
// Communicator::getGroup
// -------------------------------------------------------------
int
Communicator::getGroup(void) const
{
  return p_private->handle();
}

// -------------------------------------------------------------
// Communicator::divide
// -------------------------------------------------------------
Communicator 
Communicator::divide(int nsize) const
{
  int nprocs(size());
  int me(rank());
  // find out how many communicators need to be created
  int ngrp = nprocs/nsize;
  if (ngrp*nsize < nprocs) ngrp++;
  // evaluate how many communicators are of size nsize
  int nlarge = nprocs;
  if (ngrp*nsize > nprocs) {
    int nsmall = ngrp*nsize - nprocs;
    nlarge = nprocs - nsmall;
  }
  nlarge = nlarge/nsize;
  if (nprocs - nlarge*nsize  > nsize) nlarge++;
  // figure out what color this process is and use split functionality to
  // create new communicators
  int i;
  int color = 0;
  for (i=0; i<=me; i++) {
    if (i<nlarge*nsize) {
      if (i%nsize == 0) color++;
    } else {
      if ((i-nlarge*nsize)%(nsize-1) == 0) color++;
    }
  }
  return this->split(color);
}

// -------------------------------------------------------------
// Communicator::sync
// -------------------------------------------------------------
void 
Communicator::sync(void) const
{
  GA_Pgroup_sync(p_private->handle());
}    

} // namespace parallel
} // namespace gridpack
