// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   vector_implementation.h
 * @author William A. Perkins
 * @date   2013-05-08 08:38:39 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 26, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _vector_implementation_h_
#define _vector_implementation_h_


#include "gridpack/parallel/distributed.hpp"
#include "gridpack/utilities/uncopyable.hpp"
#include "gridpack/math/math_type.hpp"

namespace gridpack {
namespace math {

class ImplementationVisitor;

// -------------------------------------------------------------
//  class VectorImplementation
// -------------------------------------------------------------
class VectorImplementation 
  : private utility::Uncopyable,
    public parallel::Distributed
{
public:

  /// Default constructor.
  VectorImplementation(const parallel::Communicator& comm);

  /// Destructor
  ~VectorImplementation(void);

  /// Get the global length
  int size(void) const
  {
    return this->size_();
  }

  /// Get the local length
  int local_size(void) const
  {
    return this->local_size_();
  }

  // /// Set an individual element
  // void set_element(const int& i, const complex_type& x)
  // {
  //   this->set_element_(i, x);
  // }

  // /// Set an several elements
  // void set_elements(cont int& n, const int *i, const complex_type *x)
  // {
  //   this->set_elements_(n, i, x);
  // }

  // /// Add to an individual element
  // void add_element(const int& i, const complex_type& x)
  // {
  //   this->add_element_(i, x);
  // }

  // /// Add to an several elements
  // void add_elements(const int& n, const int *i, const complex_type *x)
  // {
  //   this->add_elements_(n, i, x);
  // }

  // /// Make all the elements zero
  // void zero(void)
  // {
  //   this->zero_();
  // }

  // FIXME: more ...

  /// Make this instance ready to use
  void ready(void)
  {
    this->ready_();
  }

  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor)
  {
    this->accept_(visitor);
  }

protected:

  /// Get the global vector length
  virtual int size_(void) const = 0;

  /// Get the size of the vector local part
  virtual int local_size_(void) const = 0;

  // /// Set an individual element (specialized)
  // virtual void set_element_(const int& i, const complex_type& x) = 0;

  // /// Set an several elements (specialized)
  // virtual void set_elements_(cont int& n, const int *i, const complex_type *x) = 0;

  // /// Add to an individual element (specialized)
  // virtual void add_element_(const int& i, const complex_type& x) = 0;

  // /// Add to an several elements (specialized)
  // virtual void add_elements_(const int& n, const int *i, const complex_type *x) = 0;

  // /// Make all the elements zero (specialized)
  // virtual void zero_(void) = 0;

  // FIXME: more ...
  /// Make this instance ready to use
  virtual void ready_(void) = 0;

  /// Allow visits by implementation visitors
  virtual void accept_(ImplementationVisitor& visitor) = 0;

};

} // namespace math
} // namespace gridpack



#endif
