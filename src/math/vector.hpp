// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   vector.h
 * @author William A. Perkins
 * @date   2013-05-10 12:12:59 d3g096
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
 * The values in a Vector must be set by the user.  Values are not
 * initialized.
 * 
 * When a Vector is instantiated it is ready to be filled using calls
 * to methods like ::set_element.  When the Vector is filled, all
 * processors must be notified that it is ready to use with a call to
 * ::ready().
 *
 * 
 */

class Vector 
  : public parallel::Distributed,
    public utility::Uncopyable
{
public:

  /// Default constructor.
  Vector(const parallel::Communicator& comm, const int& length);

  /// Destructor
  ~Vector(void);

  /// Get the global length
  int size(void) const
  {
    return p_vector_impl->size();
  }

  /// Get the local length
  int local_size(void) const
  {
    return p_vector_impl->local_size();
  }

  /// Get the local min/max global indexes
  void local_index_range(int& lo, int& hi) const
  {
    return p_vector_impl->local_index_range(lo, hi);
  }

  /// Set an individual element
  void set_element(const int& i, const complex_type& x)
  {
    p_vector_impl->set_element(i, x);
  }

  /// Set an several elements
  void set_elements(const int& n, const int *i, const complex_type *x)
  {
    p_vector_impl->set_elements(n, i, x);
  }

  // /// Add to an individual element
  // void add_element(const int& i, const complex_type& x)
  // {
  //   p_vector_impl->add_element(i, x);
  // }

  // /// Add to an several elements
  // void add_elements(const int& n, const int *i, const complex_type *x)
  // {
  //   p_vector_impl->add_elements(n, i, x);
  // }

  /// Get an individual element
  void get_element(const int& i, complex_type& x) const
  {
    p_vector_impl->get_element(i, x);
  }

  /// Get an several elements
  void get_elements(const int& n, const int *i, complex_type *x) const
  {
    p_vector_impl->get_elements(n, i, x);
  }

  /// Make all the elements zero
  void zero(void)
  {
    p_vector_impl->zero();
  }

  /// Make all the elements the specified value
  void fill(const complex_type& v)
  {
    p_vector_impl->fill(v);
  }

  // FIXME more ...

  /// Make this instance ready to use
  void ready(void)
  {
    p_vector_impl->ready();
  }

  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor)
  {
    p_vector_impl->accept(visitor);
  }

  void accept(ConstImplementationVisitor& visitor) const
  {
    p_vector_impl->accept(visitor);
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
  // friend Vector *reorder(const Vector& A, const Reordering& r);
  friend Vector *gridpack::math::clone(const Vector& from);


protected:
  
  /// Where stuff really happens
  boost::scoped_ptr<VectorImplementation> p_vector_impl;

  /// Constuct with an existing implementation
  explicit Vector(VectorImplementation *vimpl);
};


} // namespace utility
} // namespace gridpack

#endif
