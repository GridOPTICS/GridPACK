// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   named.h
 * @author William A. Perkins
 * @date   2013-10-09 13:44:54 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _named_h_
#define _named_h_

#include <string>

namespace gridpack {
namespace utility {

// -------------------------------------------------------------
//  class Named
// -------------------------------------------------------------
class Named {
public:

  /// Default constructor.
  Named(void);

  /// Protected copy constructor to avoid unwanted copies.
  Named(const Named& old);

  /// Destructor
  virtual ~Named(void);

  /// Get this instance's name
  std::string name(void) const
  {
    return name_;
  }

  /// Set this instance's name
  void name(const std::string& s)
  {
    name_ = s;
  }

protected:

  /// The name of this instance
  std::string name_;

};

} // namespace utility
} // namespace gridpack

#endif
