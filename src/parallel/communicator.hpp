// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   communicator.hpp
 * @author William A. Perkins
 * @date   2014-02-04 14:19:02 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------

#ifndef _communicator_hpp_
#define _communicator_hpp_

#include <boost/mpi/communicator.hpp>
#include <ga.h>

namespace gridpack {
namespace parallel {

// -------------------------------------------------------------
//  class Communicator
// -------------------------------------------------------------
class Communicator {
public:

  /// Default constructor (world)
  Communicator(void)
    : pComm()
  {
    p_GAgroup = GA_Pgroup_get_world();
  }

  /// Copy constructor.
  Communicator(const Communicator& old)
    : pComm(old.pComm)
  {
    p_GAgroup = old.p_GAgroup;
  }

  /// Construct with a Boost communicator
  /// This will only work if all processes on the world communicator
  /// call this function
  Communicator(const boost::mpi::communicator& comm)
    : pComm(comm)
  {
    int i, ierr, me, nprocs;
    me = comm.rank();
    nprocs = comm.size();
    int gaWrld = GA_Pgroup_get_world();
    int defGrp = GA_Pgroup_get_default();
    GA_Pgroup_set_default(gaWrld);
    // Construct a list of procs in the world group
    int gaSrc[nprocs];
    int gaDest[nprocs];
    for (i=0; i<nprocs; i++) {
      gaSrc[i] = 0;
      gaDest[i] = 0;
    }
    gaSrc[me] = GA_Nodeid();
    ierr = MPI_Allreduce(gaSrc,gaDest,nprocs,MPI_INT,MPI_SUM,comm);
    p_GAgroup = GA_Pgroup_create(gaDest,nprocs);
    GA_Pgroup_sync(p_GAgroup);
    GA_Pgroup_set_default(defGrp);
  }

  Communicator(const MPI_Comm& comm)
    : pComm(comm, boost::mpi::comm_duplicate)
  {
    int i, ierr, me, nprocs;
    ierr = MPI_Comm_rank(comm, &me);
    ierr = MPI_Comm_size(comm, &nprocs);
    int gaWrld = GA_Pgroup_get_world();
    int defGrp = GA_Pgroup_get_default();
    GA_Pgroup_set_default(gaWrld);
    // Construct a list of procs in the world group
    int gaSrc[nprocs];
    int gaDest[nprocs];
    for (i=0; i<nprocs; i++) {
      gaSrc[i] = 0;
      gaDest[i] = 0;
    }
    gaSrc[me] = GA_Nodeid();
    ierr = MPI_Allreduce(gaSrc,gaDest,nprocs,MPI_INT,MPI_SUM,comm);
    p_GAgroup = GA_Pgroup_create(gaDest,nprocs);
    GA_Pgroup_sync(p_GAgroup);
    GA_Pgroup_set_default(defGrp);
  }

  /// Destructor
  ~Communicator(void)
  {
    if (p_GAgroup != GA_Pgroup_get_world()){
      GA_Pgroup_destroy(p_GAgroup);
    }
  }

  /// Get the size of this communicator
  int size(void) const
  {
    return pComm.size();
  }

  /// Get this process's rank in this communicator
  int rank(void) const
  {
    return pComm.rank();
  }

  /// cast to MPI communicator
  operator MPI_Comm() const
  {
    return static_cast<MPI_Comm>(pComm);
  }

  /// cast to boost communicator
  operator boost::mpi::communicator() const
  {
    return pComm;
  }

  void barrier() const
  {
    pComm.barrier();
  }

  /// Split this instance into several
  Communicator split(int color) const
  {
    // start by creating process groups
    int i, ierr, me, nprocs;
    ierr = MPI_Comm_rank(pComm, &me);
    ierr = MPI_Comm_size(pComm, &nprocs);
    int defGrp = GA_Pgroup_get_default();
    GA_Pgroup_set_default(p_GAgroup);
    // Construct a list of procs in the world group
    int gaColor[nprocs];
    for (i=0; i<nprocs; i++) {
      gaColor[i] = 0;
    }
    gaColor[me] = color;
    GA_Pgroup_igop(p_GAgroup,gaColor,nprocs,"+");
    int ncolor = 0;
    for (i=0; i<nprocs; i++) {
      if (gaColor[i] == color) ncolor++;
    }
    int group[ncolor];
    ncolor = 0;
    for (i=0; i<nprocs; i++) {
      if (gaColor[i] == color) {
        group[ncolor] = i;
        ncolor++;
      }
    }
    int GAgroup = GA_Pgroup_create(group,ncolor);
    GA_Pgroup_set_default(defGrp);
    GA_Sync();
    Communicator ret;
    ret.pComm = pComm.split(color);
    ret.p_GAgroup = GAgroup;
    return ret;
  }

  const boost::mpi::communicator& getCommunicator(void) const
  {
    return pComm;
  }

  boost::mpi::communicator& getCommunicator(void)
  {
    return pComm;
  }

  int getGroup(void) const
  {
    return p_GAgroup;
  }

protected:
  
  /// The Boost communicator 
  boost::mpi::communicator pComm;
  int p_GAgroup;


};


} // namespace gridpack
} // namespace parallel

#endif

