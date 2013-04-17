// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   reordering_implementation.hpp
 * @author William A. Perkins
 * @date   Wed Apr 17 08:54:30 2013
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
// Created April 17, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _reordering_implementation_hpp_
#define _reordering_implementation_hpp_

// -------------------------------------------------------------
//  class ReorderingImplementation
// -------------------------------------------------------------
class ReorderingImplementation 
  : public parallel::Distributed,
    protected utility::SanityInterface,
    private utility::UnCopyable
{
public:

  /// Default constructor.
  ReorderingImplementation(const parallel::Distribution& dist, 
                           const int& current_local_size, const int& future_local_size);

  /// Destructor
  ~ReorderingImplementation(void);

  /// Make this instance ready to accept values
  void fill_begin(void) 
  {
    this->fill_begin_();
  }

  /// Indicate this instance is done accepting values (checks sanity)
  void fill_end(void) 
  {
    this->fill_end_();
  }

  /// Add a single index pair
  void add(const int& gidx_old, const int& gidx_new);

  /// Add a number of index pairs
  void add(const int& n, const int *gidx_old, const int *gidx_new);

  /// FIXME: query functions?
  
protected:

  

};



#endif
