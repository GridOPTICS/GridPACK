// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   vector.h
 * @author William A. Perkins
 * @date   Mon Mar 25 12:38:25 2013
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

#ifndef _vector_h_
#define _vector_h_

#include <boost/scoped_ptr.hpp>
#include "gridpack/parallel/distributed.hpp"
#include "gridpack/utility/uncopyable.hpp"
#include "gridpack/math/vector_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class Vector
// -------------------------------------------------------------
/// A parallel or serial vector of real values
/**
 * This class uses the Pimpl idiom for implementation in order so the
 * interface is completely free of the underlying library.  If
 * constructed with a parallel environment with only one process, a
 * serial vector is created, otherwise it's parallel. 
 * 
 */

class Vector 
  : public parallel::Distributed,
    public utility::UnCopyable
{
public:

  /// Default constructor.
  Vector(const parallel::Distribution& dist, const int& length);

  /// Destructor
  ~Vector(void);

  /// Set an individual element
  void set_element(const int& i, const double& x)
  {
    vector_impl_->set_element(i, x);
  }

  /// Set an several elements
  void set_elements(cont int& n, const int *i, const double *x)
  {
    vector_impl_->set_elements(n, i, x);
  }

  /// Add to an individual element
  void add_element(const int& i, const double& x)
  {
    vector_impl_->add_element(i, x);
  }

  /// Add to an several elements
  void add_elements(const int& n, const int *i, const double *x)
  {
    vector_impl_->add_elements(n, i, x);
  }

  /// Make all the elements zero
  void zero(void)
  {
    vector_impl_->zero();
  }

  // FIXME more ...

protected:
  
  boost::scoped_ptr<VectorImplementation> vector_impl_;
};


} // namespace utility
} // namespace gridpack

#endif
