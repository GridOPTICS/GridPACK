// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   configurable.h
 * @author William A. Perkins
 * @date   2013-10-09 13:42:39 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _configurable_h_
#define _configurable_h_

#include <string>
#include "gridpack/utility/configuration.hpp"

namespace gridpack {
namespace utility {

// -------------------------------------------------------------
//  class Configurable
// -------------------------------------------------------------
/**
 * This serves as the base for any class that has run-time
 * options. The options required by a instance of this base class are
 * identified by a string and may have any (or one of many) types.  
 *
 * Instances have two paths, one for default values for all instances
 * and one for specific instances.  The instance is created with the
 * former and it cannot be changed by any public means.  The latter
 * may be created created by prepending a name or path to the default.
 * For example, suppose a "linear solver" is a subclass, and it needs
 * a "tolerance".  The default path might be:
 * 
 *   "linear solver"/"tolerance" = 1.0e-06
 *
 * which would be used by all instances.  Say there is "nonlinear
 * solver" that uses a "linear solver" instance (both subclasses of
 * Configurable), a path such as 
 *
 *   "nonlinear solver"/"linear solver"/"tolerance" = 1.0e-03
 *
 * would distinguish it from the default.
 *
 * I think this is just a wrapper on boost::property_tree. 
 *
 */
class Configurable {
public:

  /// Copy constructor
  Configurable(const Configurable& old) : path_(old.path_) {}

  /// Destructor
  virtual ~Configurable(void) {}

  /// Get this instance's configuration path string
  std::string config_path(void) const
  {
    return path_;
  }

  /// Set this instance's configuration path string
  void config_path(const std::string& s)
  {
    path_ = s;
  }

  /// Prepend a string to this instance's path
  void config_prepend(const std::string& s);

  /// Set a particular property value
  template<typename T> void set_param(const std::string& key, const T& value);

  /// Get a particular property value
  template<typename T> T get_param(const std::string& key) const;

protected:

  /// Default constructor (only for children)
  Configurable(void) : path_() {}

  /// Construct with a specified path (only for children)
  Configurable(const std::string& p) : path_(p) {}

private:

  std::string path_;

};

} // namespace utility
} // namespace gridpack




#endif
