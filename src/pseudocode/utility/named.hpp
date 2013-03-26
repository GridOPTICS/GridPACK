// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   named.h
 * @author William A. Perkins
 * @date   Mon Mar 25 11:11:41 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 25, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

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
  ~Named(void);

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
