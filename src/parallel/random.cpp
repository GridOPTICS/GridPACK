// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   random.cpp
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

#include <ga.h>
#include <cstdlib>
#include <math.h>

#include "gridpack/parallel/random.hpp"

// -------------------------------------------------------------
//  class Random
// -------------------------------------------------------------

namespace gridpack {
namespace random {
/**
 * Default constructor
 */
Random::Random() : p_generator(42), p_uni_dist(0,1), p_uni(p_generator, p_uni_dist)
{
  p_iset = 0;
  p_gset = 0.0;
#if 0
  p_rand_max_i = 1.0/static_cast<double>(RAND_MAX);
#endif
}

/**
 * Initialize random number generator with a seed
 * @param seed random number generator initialization
 */
Random::Random(int seed) : p_generator(42), p_uni_dist(0,1), p_uni(p_generator, p_uni_dist)
{
  p_iset = 0;
  p_gset = 0.0;
  if (seed < 0) seed = -seed;
  seed += GA_Nodeid();
  p_generator.seed(static_cast<unsigned int>(seed));
#if 0
  p_rand_max_i = 1.0/static_cast<double>(RAND_MAX);
  srand(static_cast<unsigned int>(seed));
#endif
}

/**
 * Default destructor
 */
Random::~Random(void)
{
}

/**
 * Reinitialize random number generator with a new seed
 * @param seed random number generator initialization
 */
void Random::Random::seed(int seed)
{
  if (seed < 0) seed = -seed;
  seed += GA_Nodeid();
  p_generator.seed(static_cast<unsigned int>(seed));
}

/**
 * Return a double precision random number in the range [0,1]
 */
double Random::Random::drand(void)
{
  return p_uni();
#if 0
  return rand()*p_rand_max_i;
#endif
}

/**
 * Return a double precision random number from a gaussian distribution with
 * unit variance (p(x) = exp(-x*x/2)/sqrt(2*pi))
 */
double Random::Random::grand(void)
{
  if (p_iset == 0) {
    double r = 2.0;
    double x, y;
    while (r >= 1.0) {
      x = 2.0*drand() - 1.0;
      y = 2.0*drand() - 1.0;
      r = x*x + y*y;
    }
    double logr = sqrt(-2.0*log(r)/r);
    p_gset = x * logr;
    p_iset = 1;
    return y * logr;
  } else {
    p_iset = 0;
    return p_gset;
  }
  return 0.0;
}

}  // random
}  // gridpack
