// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   sanity_interface.h
 * @author William A. Perkins
 * @date   Mon Mar 25 10:56:55 2013
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
