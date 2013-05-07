// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   vector.h
 * @author William A. Perkins
 * @date   2013-05-03 10:05:08 d3g096
 * 
 * @brief  Declaration of the Vector class
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _vector_h_
#define _vector_h_

#include <boost/scoped_ptr.hpp>
#include "gridpack/parallel/distributed.hpp"
#include "gridpack/utilities/uncopyable.hpp"
#include "gridpack/math/vector_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class Vector
// -------------------------------------------------------------
/// A parallel or serial vector of values
/**
 * This class uses the Pimpl idiom for implementation  so the
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
  void set_element(const int& i, const complex_type& x)
  {
    vector_impl_->set_element(i, x);
  }

  /// Set an several elements
  void set_elements(cont int& n, const int *i, const complex_type *x)
  {
    vector_impl_->set_elements(n, i, x);
  }

  /// Add to an individual element
  void add_element(const int& i, const complex_type& x)
  {
    vector_impl_->add_element(i, x);
  }

  /// Add to an several elements
  void add_elements(const int& n, const int *i, const complex_type *x)
  {
    vector_impl_->add_elements(n, i, x);
  }

  /// Get an individual element
  void get_element(const int& i, complex_type& x) const
  {
    vector_impl_->get_element(i, x);
  }

  /// Get an several elements
  void get_elements(cont int& n, const int *i, complex_type *x) const
  {
    vector_impl_->get_elements(n, i, x);
  }

  /// Make all the elements zero
  void zero(void)
  {
    vector_impl_->zero();
  }

  // FIXME more ...

  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor)
  {
    vector_impl_->accept(visitor);
  }

  // -------------------------------------------------------------
  // In-place Vector Operation Methods (change this instance)
  // -------------------------------------------------------------
  void scale(const complex_type& x);
  void add(const Vector& x);
  void copy(const Vector& x);
  void reciprocal(void);

  // -------------------------------------------------------------
  // Vector Operations (all allocate new instances)
  // -------------------------------------------------------------
  friend Vector *add(const Vector& A, const Vector& B);
  friend Vector *reorder(const Vector& A, const Reordering& r);
  friend Vector *clone(const Vector& from);


protected:
  
  boost::scoped_ptr<VectorImplementation> vector_impl_;
};


} // namespace utility
} // namespace gridpack

#endif
