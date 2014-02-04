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
  {}

  /// Copy constructor.
  Communicator(const Communicator& old)
    : pComm(old.pComm)
  {}

  /// Construct with a Boost communicator
  Communicator(const boost::mpi::communicator& comm)
    : pComm(comm)
  {}

  Communicator(const MPI_Comm& comm)
    : pComm(comm, boost::mpi::comm_duplicate)
  {}

  /// Destructor
  ~Communicator(void) {}

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
    return Communicator(pComm.split(color));
  }

  const boost::mpi::communicator& getCommunicator(void) const
  {
    return pComm;
  }

  boost::mpi::communicator& getCommunicator(void)
  {
    return pComm;
  }

protected:
  
  /// The Boost communicator 
  boost::mpi::communicator pComm;


};


} // namespace gridpack
} // namespace parallel

#endif

