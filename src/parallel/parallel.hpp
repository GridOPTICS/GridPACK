// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   parallel.hpp
 * @author William A. Perkins
 * @date   2013-05-07 09:57:47 d3g096
 * 
 * @brief  Types and routiens used to represent the parallel environment.
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _parallel_hpp_
#define _parallel_hpp_

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

namespace gridpack {
namespace parallel {

/// A parallel environment (may need a wrapper)
typedef boost::mpi::environment Environment;

/// A communicator (may need a wrapper)
typedef boost::mpi::communicator Communicator;

} // namespace parallel
} // namespace gridpack


#endif
