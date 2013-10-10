/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   configurable.h
 * @author William A. Perkins
 * @date   2013-10-08 08:31:12 d3g096
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
#include <boost/shared_ptr.hpp>
#include "gridpack/configuration/configuration.hpp"

namespace gridpack {
namespace utility {

// -------------------------------------------------------------
//  class Configurable
// -------------------------------------------------------------
/**
 * This serves as the base for any class that has run-time options
 * that are read from the globally available Configuration
 * database. An instance of this class is assumed to have a uniquely
 * key which identifies its part of the Configuration option tree.
 *
 */
class Configurable {
public:

  /// Copy constructor
  Configurable(const Configurable& old) 
    : p_key(old.p_key),
      p_configCursor(old.p_configCursor),
      p_isConfigured(false)
  {}

  /// Destructor
  virtual ~Configurable(void) {}

  /// Get this instance's configuration key
  std::string configurationKey(void) const
  {
    return p_key;
  }

  /// Set this instance's configuration key
  void configurationKey(const std::string& s) 
  {
    p_key = s;
  }

  /// Is this instance configured?
  bool isConfigured(void) const
  {
    return p_isConfigured;
  }

  /// Initialize this instance using the specified configuration property tree
  void configure(utility::Configuration::Cursor *theprops)
  {
    if (theprops != NULL) {
      p_configCursor.reset(theprops->getCursor(this->p_key));
    }
    this->p_configure(p_configCursor.get());
    p_isConfigured = true;
  }

protected:

  /// Default constructor (only for children)
  Configurable(void) 
    : p_key("bogus"), 
      p_configCursor(),
      p_isConfigured(false)
  {}

  /// Construct with a specified path (only for children)
  Configurable(const std::string& key) 
    : p_key(key),
      p_configCursor(),
      p_isConfigured(false)
  {}

  /// Specialized way to configure from property tree
  virtual void p_configure(utility::Configuration::Cursor *props) = 0;

  /// The configuration cursor used by this instance
  /**
   * The cursor given to configure() is saved in case it needs to be
   * referred to later.  
   * 
   */
  boost::shared_ptr<utility::Configuration::Cursor> p_configCursor;

private:

  /// The configuration key
  std::string p_key;

  /// Has this instance been configured
  bool p_isConfigured;

};

} // namespace utility
} // namespace gridpack




#endif
