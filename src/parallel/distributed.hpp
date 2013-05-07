// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   distributed.hpp
 * @author William A. Perkins
 * @date   2013-05-07 10:55:28 d3g096
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
//  class Distributed
// -------------------------------------------------------------
/// Serves as a base class for parallel things
/**
 * Subclasses of this class are things that exist simultaneously on
 * multiple processes.  Instantiation and destruction occurs
 * simultaneously on all processes.
 * 
 */

class Distributed {
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
  const Communicator& communicator(void) const
  {
    return communicator_;
  }

  /// Get this processor's rank
  int processor_rank(void) const
  {
    return communicator_.rank();
  }

  /// Get the size of the parallel environment
  int processor_size(void) const
  {
    return communicator_.size();
  }
};

} // namespace parallel
} // namespace gridpack

#endif
