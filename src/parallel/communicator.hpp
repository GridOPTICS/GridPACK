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
 * @date   2014-02-13 09:25:03 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------

#ifndef _communicator_hpp_
#define _communicator_hpp_

#include <boost/mpi.hpp>
#include <boost/shared_ptr.hpp>

namespace gridpack {
namespace parallel {

// forward declaration
class CommunicatorPrivate;

// -------------------------------------------------------------
//  class Communicator
// -------------------------------------------------------------
class Communicator {
public:

  /// Default constructor (world)
  Communicator(void);

  /// Copy constructor.
  Communicator(const Communicator& old);

  /// Construct with a Boost communicator
  /** 
   * This will only work if all processes on the world communicator
   * call this function (really !?!)
   * 
   * @param comm existing Communicator
   * 
   */
  Communicator(const boost::mpi::communicator& comm);

  /// Construct from an existing @c MPI_comm
  /** 
   * 
   * 
   * @param comm 
   */
  Communicator(const MPI_Comm& comm);

  /// Make this instance the same as another
  Communicator & operator=(const Communicator & rhs);

  /// Destructor
  ~Communicator(void);

  /// Get the size of this communicator
  int size(void) const
  {
    return p_comm.size();
  }

  /// Get this process's rank in this communicator
  int rank(void) const
  {
    return p_comm.rank();
  }

  /// Get this process's rank in the world communicator
  int worldRank(void) const;

  /// cast to MPI communicator
  operator MPI_Comm() const
  {
    return static_cast<MPI_Comm>(p_comm);
  }

  /// cast to boost communicator
  operator boost::mpi::communicator() const
  {
    return p_comm;
  }

  void barrier() const
  {
    p_comm.barrier();
  }

  /// Split this instance into several communicators
  Communicator split(int color) const
  {
    return Communicator(p_comm.split(color));
  }

  /// Split this instance in several communicators
  /** 
   * Each communicator created contains at most @c nsize processes
   * 
   * @param nsize 
   * 
   * @return 
   */
  Communicator divide(int nsize) const;

  /// Get the MPI communicator for this instance (const version)
  const boost::mpi::communicator& getCommunicator(void) const
  {
    return p_comm;
  }

  /// Get a (C-style) handle to for the GA process group described by this instance
  int getGroup(void) const;


  /// Sync GA process group
  void sync(void) const;

protected:
  
  /// Swap contents with another instance
  void swap(Communicator& other);

  /// The Boost communicator 
  boost::mpi::communicator p_comm;

  /// The private (GA process group) stuff
  boost::shared_ptr<CommunicatorPrivate> p_private;
};


} // namespace gridpack
} // namespace parallel

#endif

