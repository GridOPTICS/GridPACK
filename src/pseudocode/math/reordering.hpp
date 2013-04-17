// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   reordering.hpp
 * @author William A. Perkins
 * @date   Wed Apr 17 08:38:39 2013
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April 17, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _reordering_hpp_
#define _reordering_hpp_

#include <boost/scoped_ptr.hpp>
#include "gridpack/parallel/Distributed.hpp"
#include "gridpack/math/reordering_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class Reordering
// -------------------------------------------------------------
/// Encapsulates a mapping from one index set to another
/**
 * Class are used to reorder Vector or Matrix (rows and columns)
 * instances.  It defines the both mapping between the current and
 * future @em global indices and future @em local
 * ownership.
 *
 */
class Reordering 
  : public parallel::Distributed,
    private utility::UnCopyable
{
public:

  /// Constructor 
  /** 
   * The sums current and future local sizes must be equal.
   * 
   * @param dist 
   * @param local_size local size of future ownership
   * 
   * @return 
   */
  Reordering(const parallel::Distribution& dist, 
             const int& future_local_size);

  /// Destructor
  ~Reordering(void);

  /// Make this instance ready to accept values
  void fill_begin(void)
  {
    reorder_impl_->fill_begin();
  }

  /// Indicate this instance is done accepting values (checks sanity)
  void fill_end(void)
  {
    reorder_impl_->fill_end();
  }

  /// Add a single index pair
  void add(const int& gidx_old, const int& gidx_new)
  {
    reorder_impl_->add(gidx_old, gidx_new);
  }


  /// Add a number of index pairs
  void add(const int& n, const int *gidx_old, const int *gidx_new)
  {
    reorder_impl_->add(n, gidx_old, gidx_new);
  }

  /// FIXME: query functions?
  
protected: 

  /// The actual implementation
  boost::scoped_ptr<ReorderingImplementation> reorder_impl_;

};

} // namespace math
} // namespace gridpack

#endif
