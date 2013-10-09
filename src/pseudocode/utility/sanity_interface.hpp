// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   sanity_interface.h
 * @author William A. Perkins
 * @date   2013-10-09 13:45:15 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _sanity_interface_h_
#define _sanity_interface_h_

namespace gridpack {
namespace utility {

// -------------------------------------------------------------
//  class SanityInterface
// -------------------------------------------------------------
/**
 * This class serves as the base for any class that needs to provide a
 * check of itself.
 * 
 */

class SanityInterface {
public:

  /// Default constructor.
  SanityInterface(void) {}

  /// Copy constructor
  SanityInterface(const SanityInterface& old) {}

  /// Destructor
  virtual ~SanityInterface(void) {}

  /// Is this instance OK?
  bool sane(void) const
  {
    return this->sanity_check_();
  }


protected:

  /// Specialized way to check this instance
  virtual bool sanity_check_(void) const = 0;

};

} // namespace utility
} // namespace gridpack


#endif
