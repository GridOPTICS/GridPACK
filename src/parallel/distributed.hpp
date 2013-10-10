/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   distributed.hpp
 * @author William A. Perkins
 * @date   2013-09-06 14:57:27 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------

#ifndef _distributed_hpp_
#define _distributed_hpp_

#include "gridpack/parallel/parallel.hpp"

namespace gridpack {
namespace parallel {

// -------------------------------------------------------------
//  class DistributedInterface
// -------------------------------------------------------------
/// Abstract class describing a parallel thing
/**
 * A parallel thing either contains or has access to a Communicator.
 * A parallel thing needs access to this communicator, especially to
 * determine the number processors and the local processor rank.
 * 
 */
class DistributedInterface {
protected:

public:

  /// Default constructor.
  DistributedInterface(void);

  /// Copy constructor
  DistributedInterface(const DistributedInterface& old);

  /// Destructor
  virtual ~DistributedInterface(void);

  /// Get the communicator
  virtual const Communicator& communicator(void) const = 0;

  /// Get this processor's rank
  int processor_rank(void) const
  {
    return this->communicator().rank();
  }

  /// Get the size of the parallel environment
  int processor_size(void) const
  {
    return this->communicator().size();
  }
};



// -------------------------------------------------------------
//  class Distributed
// -------------------------------------------------------------
/// Serves as a base class for parallel things
/**
 * Subclasses of this class are things that exist simultaneously on
 * multiple processes.  Instantiation and destruction occurs
 * simultaneously on all processes.  The communicator is available
 * directly to subclasses.
 * 
 */
class Distributed 
  : public DistributedInterface
{
protected:

  /// The parallel environmnet 
  Communicator communicator_;

public:

  /// Default constructor.
  explicit Distributed(const Communicator& comm);

  /// Copy constructor
  Distributed(const Distributed& old);

  /// Destructor
  virtual ~Distributed(void);

  /// Get the communicator
  const Communicator& communicator(void) const;
};

// -------------------------------------------------------------
//  class WrappedDistributed
// -------------------------------------------------------------
/// A distributed object that's just a wrapper around another distributed object
/**
 * This serves as a subclass to a class that wraps a Distributed class
 * instance.  In that case, it's desireable to have the wrapper class
 * act like a Distributed instance, but not duplicate the information
 * in the wrapped class.
 * 
 */
class WrappedDistributed 
  : public DistributedInterface 
{
protected:

  /// The actual distributed object
  Distributed *p_distributed;

  /// Set the wrapped object
  void p_set_distributed(Distributed *d) {
    p_distributed = d;
  }

  /// Default constructor.
  WrappedDistributed(Distributed *d = NULL);

public:

  /// Protected copy constructor to avoid unwanted copies.
  WrappedDistributed(const WrappedDistributed& old);

  /// Destructor
  ~WrappedDistributed(void);

  /// Get the communicator
  const Communicator& communicator(void) const;
};



} // namespace parallel
} // namespace gridpack

#endif
