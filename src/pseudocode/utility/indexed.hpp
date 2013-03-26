// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   indexed.h
 * @author William A. Perkins
 * @date   Mon Mar 25 10:48:31 2013
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
