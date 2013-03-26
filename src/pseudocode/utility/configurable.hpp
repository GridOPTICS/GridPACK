// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   configurable.h
 * @author William A. Perkins
 * @date   Mon Mar 25 11:17:21 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 25, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
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
 * options. The options required by a nistance of this base class are
 * identified by a string.
 * 
 */

class Configurable {
public:

  /// Default constructor (only for children)
  Configurable(void) : path_() {}

  /// Construct with a specified path (only for children)
  Configurable(const std::string& p) : path_(p) {}

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
protected:


private:

  std::string path_;

};

} // namespace utility
} // namespace gridpack




#endif
