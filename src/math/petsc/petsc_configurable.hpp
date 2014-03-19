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
 * @file   petsc_configurable.hpp
 * @author William A. Perkins
 * @date   2014-03-19 09:25:03 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#ifndef _petsc_configurable_hpp_
#define _petsc_configurable_hpp_

#include <gridpack/parallel/communicator.hpp>
#include <gridpack/configuration/configuration.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScConfigurable
// -------------------------------------------------------------
class PETScConfigurable 
{
public:

  /// Default constructor.
  PETScConfigurable(const parallel::Communicator& comm);

  /// Destructor
  ~PETScConfigurable(void);

  /// The common way to build a PETSc thing
  void build(utility::Configuration::CursorPtr props);

protected:

  /// Specialized way to build the PETSc parts of this instance
  virtual void p_build(const std::string& option_prefix) = 0;

private:

  /// The configuration key identifying an options prefix
  static const std::string p_prefixKey;

  /// The configuration key identifying options
  static const std::string p_optionsKey;

  /// Generate a random PETSc options prefix (for the local process)
  static std::string p_generatePrefix(void);

  /// Generate a random PETSc options prefix (for all processes)
  static std::string p_generatePrefix(const parallel::Communicator& comm);

  /// The communicator involved 
  parallel::Communicator p_comm;

  /// The option prefix used by this instance
  std::string p_prefix;

  /// The set of options added to the PETSc database
  std::string p_loadedOptions;
  
  /// Load PETSc options from configuration 
  void p_processOptions(utility::Configuration::CursorPtr props);

  /// Unload PETSc options 
  void p_unprocessOptions(void);
};


} // namespace math
} // namespace gridpack


#endif
