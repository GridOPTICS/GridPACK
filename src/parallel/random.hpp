// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   random.hpp
 * @author Bruce Palmer
 * @date   2015-03-03 09:25:03 d3g293
 * 
 * @brief  
 * This is a wrapper for a random number generator. The current implementation
 * relies on the standard random number generator in C++ and all the caveats
 * that apply to default random number generators should be noted.
 * 
 */

// -------------------------------------------------------------

#ifndef _random_hpp_
#define _random_hpp_

#include <cstdlib>

namespace gridpack {
namespace random {

// -------------------------------------------------------------
//  class GlobalIndexHashMap
// -------------------------------------------------------------
class Random {
public:

  /**
   * Default constructor
   */
  Random();

  /**
   * Initialize random number generator with a seed
   * @param seed random number generator initialization
   */
  Random(int seed);

  /**
   * Default destructor
   */
  ~Random(void);

  /**
   * Reinitialize random number generator with a new seed
   * @param seed random number generator initialization
   */
  void seed(int seed);

  /**
   * Return a double precision random number in the range [0,1]
   */
  double drand(void);

  /**
   * Return a double precision random number from a gaussian distribution with
   * unit variance
   */
  double grand(void);

private:

  int p_iset;
  double p_gset;
  double p_rand_max_i;
};


} // namespace gridpack
} // namespace random

#endif

