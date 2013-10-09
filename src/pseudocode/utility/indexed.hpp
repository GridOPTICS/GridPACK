// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   indexed.h
 * @author William A. Perkins
 * @date   2013-10-09 13:44:43 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _indexed_h_
#define _indexed_h_

namespace gridpack {
namespace utility {

// -------------------------------------------------------------
//  class Indexed
// -------------------------------------------------------------
class Indexed {
public:

  /// The index type
  /**
   * To add some flexibility. This may need to be a 64-bit integer or
   * unsigned or whatever.
   * 
   */
  typedef int IndexType;

  /// Default constructor.
  Indexed(void) 
    : local_index_(bogus_index_), global_index_(bogus_index_)
  {}

  /// Protected copy constructor to avoid unwanted copies.
  Indexed(const Indexed& old)
    : local_index_(old.local_index_),
      global_index_(old.global_index_)
  {}

  /// Destructor
  virtual ~Indexed(void) {};

  /// Get the local index of this instance
  IndexType get_local_index(void) const 
  {
    return local_index_;
  }

  /// Get the local index of this instance
  void set_local_index(const IndexType& l)
  {
    local_index_ = l;
  }

  /// Get the global index of this instance
  IndexType get_global_index(void) const 
  {
    return global_index_;
  }

  /// Set the global index of this instance
  void set_global_index(const IndexType& g)
  {
    global_index_ = g;
  }

protected:

  static const bogus_index_;
  IndexType local_index_;
  IndexType global_index_;

};

} // namespace utility
} // namespace gridpack


#endif
