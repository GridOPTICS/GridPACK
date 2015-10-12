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
 * @file   cplex_optimizer_implementation.hpp
 * @author Yilin Fang 
 * @date   2015-09-28 12:18:13 d3m045
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _cplex_optimizer_implementation_hpp_
#define _cplex_optimizer_implementation_hpp_

#include "lpfile_optimizer_implementation.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class CPlexOptimizerImplementation
// -------------------------------------------------------------
class CPlexOptimizerImplementation 
  : public LPFileOptimizerImplementation
{
public:

  /// Default constructor.
  CPlexOptimizerImplementation(const parallel::Communicator& comm);

  /// Destructor
  ~CPlexOptimizerImplementation(void);

protected:

  /// Do the problem (specialized)
  void p_solve(const p_optimizeMethod& m);

};

} // namespace optimization
} // namespace gridpack


#endif
