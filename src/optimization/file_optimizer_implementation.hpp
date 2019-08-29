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
 * @file   file_optimizer_implementation.hpp
 * @author William A. Perkins
 * @date   2016-12-08 14:59:09 d3g096
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

#ifndef _file_optimizer_implementation_hpp_
#define _file_optimizer_implementation_hpp_

#include <iosfwd>
#include <string>
#include "optimizer.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class FileOptimizerImplementation
// -------------------------------------------------------------
/// An abstracte base class for optimizers that require a temporary input file
class FileOptimizerImplementation
  : public OptimizerImplementation
{
public:

  /// Default constructor.
  FileOptimizerImplementation(const parallel::Communicator& comm);

  /// Destructor
  virtual ~FileOptimizerImplementation(void);
  
protected:

  /// Output file name
  std::string p_outputName;

  /// Try to run the file?
  bool p_runMaybe;

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props);

  /// Specialized way to set file name
  void p_setFilename(std::string file);

  /// Open a stream to a new temporary file
  virtual std::string p_temporaryFileName(void);

  /// Write an LP file to the specified stream
  virtual void p_write(const p_optimizeMethod& m, std::ostream& out) = 0;

  /// Do the problem (specialized)
  void p_solve(const p_optimizeMethod& m);
  
};

} // namespace optimization
} // namespace gridpack


#endif
