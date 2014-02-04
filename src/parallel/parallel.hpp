/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   parallel.hpp
 * @author William A. Perkins
 * @date   2014-01-31 11:25:50 d3g096
 * 
 * @brief  Types and routiens used to represent the parallel environment.
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _parallel_hpp_
#define _parallel_hpp_

#include <boost/mpi.hpp>
#include <gridpack/parallel/communicator.hpp>

namespace gridpack {
namespace parallel {

/// A parallel environment (may need a wrapper)
typedef boost::mpi::environment Environment;

} // namespace parallel
} // namespace gridpack


#endif
