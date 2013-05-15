// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   vector.h
 * @author William A. Perkins
 * @date   2013-05-10 14:27:16 d3g096
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

// forward declarations
class Vector;

/// A way to add to Vector's and make a new one
extern Vector *add(const Vector& A, const Vector& B);

/// A way to make a copy of a Vector
extern Vector *clone(const Vector& v);

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
    private utility::Uncopyable
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

  /// Add to an individual element
  void add_element(const int& i, const complex_type& x)
  {
    p_vector_impl->add_element(i, x);
  }

  /// Add to an several elements
  void add_elements(const int& n, const int *i, const complex_type *x)
  {
    p_vector_impl->add_elements(n, i, x);
  }

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

  /// Allow visits by implemetation vistor (no changes to this allowed)
  void accept(ConstImplementationVisitor& visitor) const
  {
    p_vector_impl->accept(visitor);
  }

  // -------------------------------------------------------------
  // In-place Vector Operation Methods (change this instance)
  // -------------------------------------------------------------

  /// Multiply all elements by the specified value
  void scale(const complex_type& x)
  {
    p_vector_impl->scale(x);
  }

  /// Add the specified vector
  /** 
   * This should throw an exception if the Communicator or length is
   * not the same. Local lengths can be different.
   * 
   * @param x 
   */
  void add(const Vector& x)
  {
    p_vector_impl->add(*(x.p_vector_impl));
  }

  /// Add the specified value to all elements
  void add(const complex_type& x)
  {
    p_vector_impl->add(x);
  }

  /// Copy the elements from the specified Vector
  /** 
   * This should throw an exception if the Communicator or length is
   * not the same. Local lengths can be different.
   * 
   * @param x 
   */
  void copy(const Vector& x)
  {
    p_vector_impl->add(*(x.p_vector_impl));
  }

  /// Replace all elements with their reciprocal
  void reciprocal(void)
  {
    p_vector_impl->reciprocal();
  }

  // -------------------------------------------------------------
  // Vector Operations (all allocate new instances)
  // -------------------------------------------------------------
  /// Add Vector instances
  friend Vector *add(const Vector& A, const Vector& B);

  // friend Vector *reorder(const Vector& A, const Reordering& r);

  /// Create a copy of the specified Vector
  friend Vector *clone(const Vector& from);


protected:
  
  /// Where stuff really happens
  boost::scoped_ptr<VectorImplementation> p_vector_impl;

  /// Constuct with an existing implementation
  explicit Vector(VectorImplementation *vimpl);
};


} // namespace utility
} // namespace gridpack

#endif
