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
 * @file   lpfile_optimizer_implementation.hpp
 * @author William A. Perkins
 * @date   2015-10-13 13:11:36 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#ifndef _lpfile_optimizer_implementation_hpp_
#define _lpfile_optimizer_implementation_hpp_

#include <iosfwd>
#include <string>
#include <map>
#include "optimizer.hpp"

namespace gridpack {
namespace optimization {



// -------------------------------------------------------------
//  class LPFileOptimizerImplementation
// -------------------------------------------------------------
class LPFileOptimizerImplementation 
  : public OptimizerImplementation
{
public:

  /// Default constructor.
  LPFileOptimizerImplementation(const parallel::Communicator& comm)
    : OptimizerImplementation(comm)
  {}

  /// Destructor
  ~LPFileOptimizerImplementation(void)
  {}

protected:

  /// Open a stream to a new temporary file
  std::string p_temporaryFileName(void);

  /// Write an LP file to the specified stream
  virtual void p_write(const p_optimizeMethod& m, std::ostream& out);

  /// Do the problem (specialized)
  void p_solve(const p_optimizeMethod& m);
  
};


} // namespace optimization
} // namespace gridpack


#endif
