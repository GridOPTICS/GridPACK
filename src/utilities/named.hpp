//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   named.h
 * @author William A. Perkins
 * @date   2017-02-10 07:31:37 d3g096
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

#include <boost/serialization/string.hpp>
#include <boost/serialization/export.hpp>

namespace gridpack {
namespace utility {

// -------------------------------------------------------------
//  class Named
// -------------------------------------------------------------
class Named {
public:

  /// Default constructor.
  Named(void)
    : p_name()
  {}

  /// Construct with a value
  Named(const std::string& s)
    : p_name(s)
  {}

  /// Construct with a C-string
  Named(const char *s)
    : p_name(s)
  {}

  /// Protected copy constructor to avoid unwanted copies.
  Named(const Named& old)
    : p_name(old.p_name)
  {}

  /// Destructor
  virtual ~Named(void) {}

  /// Get this instance's name
  virtual std::string name(void) const
  {
    return p_name;
  }

  /// Set this instance's name
  void name(const std::string& s)
  {
    p_name = s;
  }

protected:

  /// The name of this instance
  std::string p_name;

private:

  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & p_name;
  }
};

} // namespace utility
} // namespace gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::utility::Named)

#endif
