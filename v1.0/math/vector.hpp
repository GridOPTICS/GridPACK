// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   vector.h
 * @author William A. Perkins
 * @date   2013-10-28 13:44:47 d3g096
 * 
 * @brief  Declaration of the Vector class
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _vector_h_
#define _vector_h_

#include <boost/scoped_ptr.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/math/vector_implementation.hpp>

namespace gridpack {
namespace math {

// forward declarations
class Vector;

// -------------------------------------------------------------
//  class Vector
// -------------------------------------------------------------
/// A parallel or serial vector of values
/**
 * This class encapsulates a vector of values.  
 * 
 * When a Vector is instantiated it is ready to be filled using calls
 * to methods like setElement().  When the Vector is filled, all
 * processors must be notified that it is ready to use with a call to
 * ready().  
 *
 * For efficiency, setElement(), and related calls, should only set
 * locally owned elements, even though off processor elements can be
 * set with those methods.  getElement(), and related calls, can only
 * access locally owned elements.  Off processor values can be obtained
 * by calling get_all_elements() (currently the only way).
 *
 * Element indexing is always global and 0-based.
 *
 * This class uses the Pimpl idiom for \ref VectorImplementation
 * "implementation" so the interface is completely free of the
 * underlying library.  This class simply provides an interface to a
 * specific \ref VectorImplementation "implementation".
 * 
 */
class Vector 
  : public parallel::WrappedDistributed,
    private utility::Uncopyable
{
public:

  /// Default constructor.
  /** 
   * @e Collective.
   *
   * A vector must be instantiated simulutaneously on all processes
   * involved in the specified \ref parallel::Communicator
   * "communicator". Each process in the communicator will own the
   * number of elements requested.
   * 
   * @param comm parallel environment 
   * @param local_length number of elements to be assigned to this process
   * 
   * @return empty vector instance
   */
  Vector(const parallel::Communicator& comm, const int& local_length);

  /// Destructor
  /** 
   * @e Collective.
   *
   * A vector must be destroyed simulutaneously on all processes
   * involved in the \ref parallel::Communicator "communicator" used
   * to instantiate it.
   */
  ~Vector(void);

  /// Get the global vector length
  /** 
   * @e Local.
   * 
   * 
   * 
   * 
   * @return global vector length
   */
  int size(void) const
  {
    return p_vector_impl->size();
  }

  /// Get the number of locally owned vector elements
  /** 
   * @e Local.
   * 
   * 
   * 
   * 
   * @return local vector length
   */
  int localSize(void) const
  {
    return p_vector_impl->localSize();
  }

  /** 
   * @deprecated does not meet coding standards, use localSize()
   * 
   * 
   * @return 
   */
  int local_size(void) const
  {
    return this->localSize();
  }
  

  /// Get the min/max global indexes of locally owned elements
  /** 
   * @e Local.
   * 
   * The minimum index in the first global (0-based) index owned by
   * the local process. The maximum is one more than the last global
   * index owned. An example of usage:
   * 
     \code{.cpp}
     int lo, hi;
     my_vector.local_index_range(lo, hi);
     for (int i = lo; i < hi; ++i) {
       ComplexType x;
       x = ...;
       v.setElement(i, x);
     }
     \endcode
   * 
   * @param lo first (0-based) index of locally owned elements
   * @param hi one more than the last (0-based) index of locally owned elements
   */
  void localIndexRange(int& lo, int& hi) const
  {
    return p_vector_impl->localIndexRange(lo, hi);
  }

  /** 
   * @deprecated does not meet coding standards, use localIndexRange()
   * 
   * @param lo 
   * @param hi 
   */
  void local_index_range(int& lo, int& hi) const
  {
    return this->localIndexRange(lo, hi);
  }


  /// Set an individual element
  /** 
   * @e Local.
   * 
   * This overwrites the value at the specified index.  ready() must
   * be called after all setElement() calls and before using the
   * vector.
   * 
   * @param i element global (0-based) index 
   * @param x value to place in vector
   */
  void setElement(const int& i, const ComplexType& x)
  {
    p_vector_impl->setElement(i, x);
  }

  /** 
   * @deprecated does not meet coding standards, use setElement
   * 
   * @param i 
   * @param x 
   */
  void set_element(const int& i, const ComplexType& x)
  {
    this->setElement(i, x);
  }

  /// Set several elements
  /** 
   * @e Local.
   * 
   * This places (overwrites) several elements, with arbitrary
   * indexes, in the vector.  ready() must be called after all
   * setElement() calls and before using the vector.
   * 
   * @param n number of elements to place in vector
   * @param i pointer to an array of @c n global (0-based) indexes
   * @param x pointer to an arry 
   */
  void setElements(const int& n, const int *i, const ComplexType *x)
  {
    p_vector_impl->setElements(n, i, x);
  }

  /** 
   * @deprecated does not meet coding standards, use setElements()
   * 
   * @param n 
   * @param i 
   * @param x 
   */
  void set_elements(const int& n, const int *i, const ComplexType *x)
  {
    this->setElements(n, i, x);
  }

  /// Set a range of elements (lo to hi-1)
  /** 
   * @e Local.
   * 
   * An example that fills the locally owned part of the vector:
   * \code{.cpp}
   * Vector v(...);
   * int lo, hi;
   * v.local_index_range(lo, hi);
   * std::vector<ComplexType> x(v.local_size())
   * // fill x with appropriate values
   * v.setElement_range(lo, hi, &x[0]);
   * \endcode
   * 
   * @param lo lowest global (0-based) index to fill
   * @param hi one more than the highest global (0-based) index to fill
   * @param x array of hi - lo values
   */
  void setElementRange(const int& lo, const int& hi, ComplexType *x)
  {
    p_vector_impl->setElementRange(lo, hi, x);
  }

  /** 
   * @deprecated does not meet coding standards, use setElementRange()
   * 
   * @param lo 
   * @param hi 
   * @param x 
   */
  void set_element_range(const int& lo, const int& hi, ComplexType *x)
  {
    this->setElementRange(lo, hi, x);
  }

  /// Add to an individual element
  /** 
   * @e Local.
   * 
   * 
   * 
   * @param i 
   * @param x 
   */
  void addElement(const int& i, const ComplexType& x)
  {
    p_vector_impl->addElement(i, x);
  }

  /** 
   * @deprecated does not meet coding standards, use addElement()
   * 
   * @param i 
   * @param x 
   */
  void add_element(const int& i, const ComplexType& x)
  {
    this->addElement(i, x);
  }

  /// Add to an several elements
  /** 
   * @e Local.
   * 
   * 
   * 
   * @param n 
   * @param i 
   * @param x 
   */
  void addElements(const int& n, const int *i, const ComplexType *x)
  {
    p_vector_impl->addElements(n, i, x);
  }

  /** 
   * @deprecated does not meet coding standards, use addElements()
   * 
   * @param n 
   * @param i 
   * @param x 
   */
  void add_elements(const int& n, const int *i, const ComplexType *x)
  {
    this->addElements(n, i, x);
  }

  /// Get an individual element
  /** 
   * @e Local.
   * 
   * Only local values may be retreived with this method.  
   * 
   * @param i element's global, 0-based index
   * @param x element's value
   */
  void getElement(const int& i, ComplexType& x) const
  {
    p_vector_impl->getElement(i, x);
  }

  /** 
   * @deprecated does not meet coding standards, use getElement()
   * 
   * @param i 
   * @param x 
   */
  void get_element(const int& i, ComplexType& x) const
  {
    this->getElement(i, x);
  }

  /// Get an several elements
  /** 
   * @e Local.
   * 
   * Only local values may be retreived with this method.  
   * 
   * 
   * @param n number of elements to get
   * @param i array of @c n global, 0-based indexes
   * @param x array of @c n values
   */
  void getElements(const int& n, const int *i, ComplexType *x) const
  {
    p_vector_impl->getElements(n, i, x);
  }

  /** 
   * @deprecated does not meet coding standards, use getElements()
   * 
   * @param n 
   * @param i 
   * @param x 
   */
  void get_elements(const int& n, const int *i, ComplexType *x) const
  {
    this->getElements(n, i, x);
  }


  /// Get a range of elements
  /** 
   * @e Local.
   *
   * The elements whose indexes range from @c lo to @c hi-1 are
   * retreived.  Only local values may be retreived with this method.
   * The length of the @c x array must be at least @c "hi - lo"
   * 
   * @param lo lowest element index
   * @param hi one more than the highest element index
   * @param x existing array of values to be filled with element values
   */
  void getElementRange(const int& lo, const int& hi, ComplexType *x) const
  {
    p_vector_impl->getElementRange(lo, hi, x);
  }

  /** 
   * @deprecated does not meet coding standards, use getElementRange()
   * 
   * @param lo 
   * @param hi 
   * @param x 
   */
  void get_element_range(const int& lo, const int& hi, ComplexType *x) const
  {
    this->getElementRange(lo, hi, x);
  }

  /// Get all of the vector elements
  /** 
   * @e Collective.
   *
   * This is an all gather operation and consequently will be slow.
   * All of this vector's elements are gathered on all processes and
   * placed in the in the specified array.
   * 
   * 
   * @param x array of size() length to be filled with values
   */
  void getAllElements(ComplexType *x) const
  {
    p_vector_impl->getAllElements(x);
  }

  /** 
   * @deprecated does not meet coding standards, use getAllElements()
   * 
   * @param x 
   */
  void get_all_elements(ComplexType *x) const
  {
    this->getAllElements(x);
  }

  /// Make all the elements zero
  /** 
   * @e Collective.
   *
   * 
   * 
   */
  void zero(void)
  {
    p_vector_impl->zero();
  }

  /// Make all the elements the specified value
  /** 
   * @e Collective.
   *
   * 
   * 
   * @param v 
   */
  void fill(const ComplexType& v)
  {
    p_vector_impl->fill(v);
  }

  /// Compute the vector L<sup>1</sup> norm
  /** 
   * @e Collective.
   *
   * The vector L<sup>1</sup> norm is computed as
   * \f[
   *   \left\| \mathbf{x} \right\|_{1} ~ = ~ \sum_{i} \left| x_{i} \right|
   * \f]
   * 
   * 
   * @return 
   */
  ComplexType norm1(void) const
  {
    return p_vector_impl->norm1();
  }

  /// Compute the vector L<sup>2</sup> norm
  /** 
   * @e Collective.
   *
   * The vector L<sup>2</sup>, or Euclidian, norm is computed as
   * \f[
   *   \left\| \mathbf{x} \right\| ~ = ~ \sqrt{\sum_{i} x_{i}^{2}}
   * \f]
   * 
   * @return 
   */
  ComplexType norm2(void) const
  {
    return p_vector_impl->norm2();
  }

  // FIXME more ...

  /// Make this instance ready to use
  /** 
   * @e Collective.
   *
   * This is used to indicate that the vector is ready to use.  This
   * must be called after @e all setElement() or add element() calls
   * and before the vector is used for any operation.
   * 
   */
  void ready(void)
  {
    p_vector_impl->ready();
  }

  //! @cond DEVDOC

  /// Allow visits by implemetation visitor
  /** 
   * 
   * 
   * @param visitor 
   */
  void accept(ImplementationVisitor& visitor)
  {
    p_vector_impl->accept(visitor);
  }

  /// Allow visits by implemetation vistor (no changes to this allowed)
  /** 
   * 
   * 
   * @param visitor 
   */
  void accept(ConstImplementationVisitor& visitor) const
  {
    p_vector_impl->accept(visitor);
  }

  //! @endcond

  /// Make an exact replica of this instance
  /** 
   * @e Collective.
   *
   * @return 
   */
  Vector *clone(void) const
  {
    VectorImplementation *pimpl_clone = p_vector_impl->clone();
    Vector *result = new Vector(pimpl_clone);
    return result;
  }

  /// Print to named file or standard output
  /** 
   * @e Collective.
   *
   * 
   *
   * The format is dependent on the specific vector implementation.
   * 
   * @param filename optional file
   */
  void print(const char* filename = NULL) const;

  /// Save, in MatLAB format, to named file (collective)
  /** 
   * @e Collective.
   *
   * 
   * 
   * @param filename 
   */
  void save(const char *filename) const;

  // -------------------------------------------------------------
  // In-place Vector Operation Methods (change this instance)
  // -------------------------------------------------------------

  /// Multiply all elements by the specified value
  /** 
   * @e Collective.
   *
   * 
   * 
   * @param x 
   */
  void scale(const ComplexType& x);

  /// Add the specified vector
  /** 
   * @e Collective.
   *
   * 
   * 
   * @param x 
   * @param scale 
   */
  void add(const Vector& x, const ComplexType& scale = 1.0);

  /// Add the specified value to all elements
  /** 
   * @e Collective.
   *
   * 
   * 
   * @param x 
   */
  void add(const ComplexType& x);

  /// Copy the elements from the specified Vector
  /** 
   * @e Collective.
   *
   * 
   * 
   * @param x 
   */
  void equate(const Vector& x);

  /// Replace all elements with their reciprocal
  /** 
   * @e Collective.
   *
   * 
   * 
   */
  void reciprocal(void);

  /// Replace all elements with its absolute value (complex magnitude) 
  void abs(void)
  {
    p_vector_impl->abs();
  }

  /// Replace all elements with their real part
  void real(void)
  {
    p_vector_impl->real();
  }

  /// Replace all elements with their imaginary part
  void imaginary(void)
  {
    p_vector_impl->imaginary();
  }

  /// Replace all elements with their complex conjugate
  void conjugate(void)
  {
    p_vector_impl->conjugate();
  }

  // friend Vector *reorder(const Vector& A, const Reordering& r);

protected:
  
  /// Where stuff really happens
  boost::scoped_ptr<VectorImplementation> p_vector_impl;

  /// Constuct with an existing implementation
  /** 
   * This constructor is used to create a new vector instance using an
   * existing VectorImplementation instance.  This is needed by
   * certain vector operations.
   * 
   * @param vimpl pointer to a new, allocated VectorImplementation instance
   * 
   * @return new vector instance using the specified implementation
   */
  explicit Vector(VectorImplementation *vimpl);

  /// Is this Vector compatible with this one, throw if not
  /** 
   * @e Local.
   *
   * This checks to see if @c x uses the same \ref
   * parallel::Communicator "communicator" and has the same global
   * length as this instance. If @c x is not compatible, an exception
   * is thrown.  @c x may have a different parallel distribution than
   * this instance.
   * 
   * @param x vector to check
   */
  void p_checkCompatible(const Vector& x) const;
};

// -------------------------------------------------------------
// Vector Operations (all allocate new instances)
// -------------------------------------------------------------

/// Add two Vector instances and put the result in a new one
/** 
 * @e Collective.
 *
 * @c A and @c B must have the global size and use the same \ref
 * parallel::Communicator "communicator", but may have different local
 * sizes.
 *
 * The resultant vector instance will have the same parallel
 * distribution as @c A.
 * 
 * @param A 
 * @param B 
 * 
 * @return pointer to new allocated Vector instance
 */
extern Vector *add(const Vector& A, const Vector& B);

/// Subtract two Vector instances and put the result in a new one
/** 
 * @e Collective.
 *
 * @c A and @c B must have the global size and use the same \ref
 * parallel::Communicator "communicator", but may have different local
 * sizes.
 *
 * The resultant vector instance will have the same parallel
 * distribution as @c A.
 * 
 * @param A 
 * @param B 
 * 
 * @return pointer to new allocated Vector instance
 */
extern Vector *subtract(const Vector& A, const Vector& B);

/// Create a vector containing the absolute value/magnitude of the specified vector
extern Vector *abs(const Vector& x);

/// Create a vector containing the real part of the specified vector
/** 
 * 
 * 
 * @param x existing vector with complex values
 *  
 * @return pointer to new vector instance containing real part of @c x
 */
extern Vector *real(const Vector& x);

/// Create a vector containing the imaginar part of the specified vector
/** 
 * 
 * 
 * @param x existing vector with complex values
 * 
 * @return pointer to new vector instance containing imaginary part of @c x
 */
extern Vector *imaginary(const Vector& x);

/// Create a vector containing the complex conjugate of @c x
/** 
 * 
 * 
 * @param x existing vector with complex values
 * 
 * @return pointer to new vector instance containing the complex conjugate of @c x
 */
extern Vector* conjugate(const Vector& x);



// -------------------------------------------------------------
// Vector Operations (results into existing instances)
// -------------------------------------------------------------

/// Add two Vector instances and put the result in an existing Vector
/** 
 * @e Collective.
 *
 * @c A, @c B, and @c result need to have the same global size and use the
 * same \ref parallel::Communicator "communicator", but may have
 * different local sizes.
 * 
 * @param A 
 * @param B 
 * @param result vector in which to place the sum of @c A and @c B
 */
extern void add(const Vector& A, const Vector& B, Vector& result);



} // namespace utility
} // namespace gridpack

#endif
