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
 * @date   2013-11-08 09:00:49 d3g096
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
//  class ConfigurableInterface
// -------------------------------------------------------------
class ConfigurableInterface {
public:

  /// Copy constructor
  ConfigurableInterface(const ConfigurableInterface& old) {};

  /// Destructor
  virtual ~ConfigurableInterface(void) {};

  /// Get this instance's configuration key
  virtual std::string configurationKey(void) const = 0;

  /// Set this instance's configuration key
  virtual void configurationKey(const std::string& s) = 0;

  /// Is this instance configured?
  virtual bool isConfigured(void) const = 0;

  /// Initialize this instance using the specified configuration property tree
  virtual void configure(utility::Configuration::CursorPtr theprops) = 0;

protected:

  /// Default constructor (only for children)
  ConfigurableInterface(void) {};

};


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
class Configurable 
  : public ConfigurableInterface
{
public:

  /// Copy constructor
  Configurable(const Configurable& old) 
    : p_key(old.p_key),
      p_configCursor(old.p_configCursor),
      p_isConfigured(false)
  {}

  /// Destructor
  ~Configurable(void) {}

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
  void configure(utility::Configuration::CursorPtr theprops)
  {
    if (theprops) {
      p_configCursor = theprops->getCursor(this->p_key);
    }
    this->p_configure(p_configCursor);
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
  virtual void p_configure(utility::Configuration::CursorPtr props) = 0;

  /// The configuration cursor used by this instance
  /**
   * The cursor given to configure() is saved in case it needs to be
   * referred to later.  
   * 
   */
  utility::Configuration::CursorPtr p_configCursor;

private:

  /// The configuration key
  std::string p_key;

  /// Has this instance been configured
  bool p_isConfigured;

};

// -------------------------------------------------------------
//  class WrappedConfigurable
// -------------------------------------------------------------
class WrappedConfigurable 
  : public ConfigurableInterface
{
public:

  /// Copy constructor
  WrappedConfigurable(const WrappedConfigurable& old) 
    : ConfigurableInterface(old), 
      p_configurable(old.p_configurable)
  {}

  /// Destructor
  ~WrappedConfigurable(void) 
  {}

  /// Get this instance's configuration key
  std::string configurationKey(void) const
  {
    return p_configurable->configurationKey();
  }

  /// Set this instance's configuration key
  void configurationKey(const std::string& s) 
  {
    p_configurable->configurationKey(s);
  }

  /// Is this instance configured?
  bool isConfigured(void) const
  {
    return p_configurable->isConfigured();
  }

  /// Initialize this instance using the specified configuration property tree
  void configure(utility::Configuration::CursorPtr theprops)
  {
    p_configurable->configure(theprops);
  }

protected:

  /// Default constructor (only for childred)
  WrappedConfigurable(void)
    : ConfigurableInterface(), p_configurable()
  {}

  /// Constructor with existing Configurable instance (only for childred)
  WrappedConfigurable(Configurable *c) 
    : ConfigurableInterface(), p_configurable(c)
  {}

  /// The wrapped Configurable instance
  Configurable *p_configurable;

  /// Set the wrapped Configurable instance
  void p_setConfigurable(Configurable *c)
  {
    p_configurable = c;
  }

};



} // namespace utility
} // namespace gridpack




#endif
