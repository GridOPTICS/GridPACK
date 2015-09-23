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
 * @file   glpk_optimizer_implementation.hpp
 * @author William A. Perkins
 * @date   2015-09-16 10:18:13 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _glpk_optimizer_implementation_hpp_
#define _glpk_optimizer_implementation_hpp_

#include "lpfile_optimizer_implementation.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class GLPKOptimizerImplementation
// -------------------------------------------------------------
class GLPKOptimizerImplementation 
  : public LPFileOptimizerImplementation
{
public:

  /// Default constructor.
  GLPKOptimizerImplementation(const parallel::Communicator& comm);

  /// Destructor
  ~GLPKOptimizerImplementation(void);

protected:

  /// Do the problem (specialized)
  void p_solve(const p_optimizeMethod& m);

};

} // namespace optimization
} // namespace gridpack


#endif
