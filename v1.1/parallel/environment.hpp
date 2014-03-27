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
 * @file   environment.hpp
 * @author William A. Perkins
 * @date   2014-02-10 14:06:28 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _environment_hpp_
#define _environment_hpp_

#include <boost/mpi.hpp>
#include "gridpack/utilities/uncopyable.hpp"

namespace gridpack {
namespace parallel {


// -------------------------------------------------------------
//  class Environment
// -------------------------------------------------------------
class Environment 
  : private utility::Uncopyable
{
public:

  /// Default constructor.
  /** 
   * 
   * 
   * @param argc number of program command line arguments
   * @param argv command line arguments
   * @param ma_stack stack size, bytes, to be allocated for Global Arrays
   * @param ma_heap heap size, bytes, to be allocated for Global Arrays
   * 
   * @return 
   */
  Environment(int& argc, char **argv,
              const long int& ma_stack = 200000,
              const long int& ma_heap = 200000);

  /// Destructor
  ~Environment(void);

protected:

  /// The (Boost) MPI environment
  boost::mpi::environment p_boostEnv;

  // FIXME: there are some static boost::mpi methods that should be
  // wrapped here

};



} // namespace parallel
} // namespace gridpack

#endif
