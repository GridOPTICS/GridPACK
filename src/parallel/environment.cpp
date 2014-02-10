// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   environment.cpp
 * @author William A. Perkins
 * @date   2014-02-10 14:07:08 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <ga++.h>
#include "environment.hpp"

namespace gridpack {
namespace parallel {

// -------------------------------------------------------------
//  class Environment
// -------------------------------------------------------------

// -------------------------------------------------------------
// Environment:: constructors / destructor
// -------------------------------------------------------------
Environment::Environment(int& argc, char **argv,
                         const long int& ma_stack, 
                         const long int& ma_heap)
  : p_boostEnv(argc, argv)
{
  GA_Initialize();
  MA_init(MT_C_CHAR, ma_stack, ma_heap);
}

Environment::~Environment(void)
{
  GA_Terminate();
}


} // namespace parallel
} // namespace gridpack

