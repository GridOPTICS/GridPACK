// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   configurable.h
 * @author William A. Perkins
 * @date   2013-10-01 11:45:34 d3g096
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
  Configurable(const Configurable& old) : p_key(old.p_key) {}

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

  /// Initialize this instance using the specified configuration property tree
  void configure(utility::Configuration::Cursor *theprops)
  {
    utility::Configuration::Cursor *props(NULL);
    if (theprops != NULL) {
      props = theprops->getCursor(this->p_key);
    }
    this->p_configure(props);
  }

protected:

  /// Default constructor (only for children)
  Configurable(void) : p_key("bogus") {}

  /// Construct with a specified path (only for children)
  Configurable(const std::string& key) : p_key(key) {}

  /// Specialized way to configure from property tree
  virtual void p_configure(utility::Configuration::Cursor *props) = 0;

private:

  /// The configuration key
  std::string p_key;

};

} // namespace utility
} // namespace gridpack




#endif
