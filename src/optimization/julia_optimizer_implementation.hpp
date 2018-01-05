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
 * @file   julia_optimizer_implementation.hpp
 * @author William A. Perkins
 * @date   2016-12-08 14:56:41 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December  7, 2016 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _julia_optimizer_implementation_hpp_
#define _julia_optimizer_implementation_hpp_

#include "file_optimizer_implementation.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class JuliaOptimizerImplementation
// -------------------------------------------------------------
/// Base class for optimizers that take the Julia language as input
class JuliaOptimizerImplementation
  : public FileOptimizerImplementation
{
public:

  /// Default constructor.
  JuliaOptimizerImplementation(const parallel::Communicator& comm);

  /// Destructor
  ~JuliaOptimizerImplementation(void);

protected:

  /// Some extra stuff to put at the beginning of the Julia script
  std::string p_preamble;

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props);

  /// Specialized way to set file name
  void p_setFilename(std::string file);

  /// Write an Julia file to the specified stream
  void p_write(const p_optimizeMethod& m, std::ostream& out);

};



} // namespace optimization
} // namespace gridpack

#endif
