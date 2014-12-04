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
 * @date   2014-10-30 14:20:46 d3g096
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
    private utility::Uncopyable,
    public BaseVectorInterface<ComplexType>
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
  explicit Vector(VectorImplementation<ComplexType> *vimpl);

  /// Destructor
  /** 
   * @e Collective.
   *
   * A vector must be destroyed simulutaneously on all processes
   * involved in the \ref parallel::Communicator "communicator" used
   * to instantiate it.
   */
  ~Vector(void);

  //! @cond DEVDOC

  //! @endcond

  /// Make an exact replica of this instance
  /** 
   * @e Collective.
   *
   * @return 
   */
  Vector *clone(void) const
  {
    VectorImplementation<ComplexType> *pimpl_clone = p_vector_impl->clone();
    Vector *result = new Vector(pimpl_clone);
    return result;
  }

  // -------------------------------------------------------------
  // In-place Vector Operation Methods (change this instance)
  // -------------------------------------------------------------

  /// Add the specified vector
  /** 
   * @e Collective.
   *
   * 
   * 
   * @param x 
   * @param scale 
   */
  void add(const Vector& x, const TheType& scale = 1.0);

  /// Add the specified value to all elements
  /** 
   * @e Collective.
   *
   * 
   * 
   * @param x 
   */
  void add(const TheType& x);

  /// Copy the elements from the specified Vector
  /** 
   * @e Collective.
   *
   * 
   * 
   * @param x 
   */
  void equate(const Vector& x);

  /// Element-by-element multiply by another Vector
  void elementMultiply(const Vector& x);

  /// ELement-by-element divide by another Vector
  void elementDivide(const Vector& x);

  // friend Vector *reorder(const Vector& A, const Reordering& r);

protected:
  
  /// Where stuff really happens
  boost::scoped_ptr< VectorImplementation<ComplexType> > p_vector_impl;

  /// Get the global vector length (specialized)
  IdxType p_size(void) const
  { p_vector_impl->size(); }

  /// Get the size of the vector local part (specialized)
  IdxType p_localSize(void) const
  { p_vector_impl->localSize(); }

  /// Get the local min/max global indexes (specialized)
  void p_localIndexRange(IdxType& lo, IdxType& hi) const
  { p_vector_impl->localIndexRange(lo, hi); }

  /// Set an individual element (specialized)
  void p_setElement(const IdxType& i, const TheType& x)
  { p_vector_impl->setElement(i, x); }

  /// Set an several elements (specialized)
  void p_setElements(const IdxType& n, const IdxType *i, const TheType *x)
  { p_vector_impl->setElements(n, i, x); }

  /// Get a range of elements (lo to hi-1) (specialized)
  void p_setElementRange(const IdxType& lo, const IdxType& hi, TheType *x)
  { p_vector_impl->setElementRange(lo, hi, x); }

  /// Add to an individual element (specialized)
  void p_addElement(const IdxType& i, const TheType& x)
  { p_vector_impl->addElement(i, x); }

  /// Add to an several elements (specialized)
  void p_addElements(const IdxType& n, const IdxType *i, const TheType *x)
  { p_vector_impl->addElements(n, i, x); }

  /// Get an individual element (specialized)
  void p_getElement(const IdxType& i, TheType& x) const
  { p_vector_impl->getElement(i, x); }

  /// Get an several elements (specialized)
  void p_getElements(const IdxType& n, const IdxType *i, TheType *x) const
  { p_vector_impl->getElements(n, i, x); }

  /// Get a range of elements (lo to hi-1) (specialized)
  void p_getElementRange(const IdxType& lo, const IdxType& hi, TheType *x) const
  { p_vector_impl->getElementRange(lo, hi, x); }

  /// Get all of vector elements (on all processes)
  void p_getAllElements(TheType *x) const
  { p_vector_impl->getAllElements(x); }

  /// Make all the elements zero (specialized)
  void p_zero(void)
  { p_vector_impl->zero(); }

  /// Fill all the elements with the specified value (specialized)
  void p_fill(const TheType& v)
  { p_vector_impl->fill(v); }

  /// Scale all elements by a single value (specialized)
  void p_scale(const TheType& x)
  { p_vector_impl->scale(x); }

  /// Compute the vector L1 norm (sum of absolute value) (specialized)
  double p_norm1(void) const
  { return p_vector_impl->norm1(); }

  /// Compute the vector L2 norm (root of sum of squares) (specialized)
  double p_norm2(void) const
  { return p_vector_impl->norm2(); }

  /// Compute the vector infinity (or maximum) norm (specialized)
  double p_normInfinity(void) const
  { return p_vector_impl->normInfinity(); }

  /// Replace all elements with its absolute value (specialized) 
  void p_abs(void)
  { p_vector_impl->abs(); }

   /// Replace all elements with their real part (specialized)
  void p_real(void)
  { p_vector_impl->real(); }
 
   /// Replace all elements with their imaginary part (specialized)
  void p_imaginary(void)
  { p_vector_impl->imaginary(); }
 
   /// Replace all elements with their complex conjugate
  void p_conjugate(void)
  { p_vector_impl->conjugate(); }
 
  /// Replace all elements with its exponential (specialized)
  void p_exp(void)
  { p_vector_impl->exp(); }

  /// Replace all elements with its reciprocal (specialized)
  void p_reciprocal(void)
  { p_vector_impl->reciprocal(); }

  /// Make this instance ready to use (specialized)
  void p_ready(void)
  { p_vector_impl->ready(); }

  /// Print to named file or standard output (specialized)
  void p_print(const char* filename = NULL) const 
  { p_vector_impl->print(filename); }

  /// Save, in MatLAB format, to named file (collective) (specialized)
  void p_save(const char *filename) const 
  { p_vector_impl->save(filename); }

  /// Load from a named file of whatever binary format the math library uses (specialized)
  void p_loadBinary(const char *filename) 
  { p_vector_impl->loadBinary(filename); }

  /// Save to named file in whatever binary format the math library uses (specialized)
  void p_saveBinary(const char *filename) const 
  { p_vector_impl->saveBinary(filename); }

  /// Allow visits by implemetation visitor
  void p_accept(ImplementationVisitor& visitor)
  {
    p_vector_impl->accept(visitor);
  }

  /// Allow visits by implemetation vistor (no changes to this allowed)
  void p_accept(ConstImplementationVisitor& visitor) const
  {
    p_vector_impl->accept(visitor);
  }

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
inline Vector *add(const Vector& A, const Vector& B)
{
  Vector *result(A.clone());
  result->add(B);
  return result;
}

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
// inline Vector *subtract(const Vector& A, const Vector& B);

/// Create a vector containing the absolute value/magnitude of the specified vector
inline Vector *abs(const Vector& x)
{
  Vector *result(x.clone());
  result->abs();
  return result;
}

/// Create a vector containing the real part of the specified vector
/** 
 * 
 * 
 * @param x existing vector with complex values
 *  
 * @return pointer to new vector instance containing real part of @c x
 */
inline Vector *real(const Vector& x)
{
  Vector *result(x.clone());
  result->real();
  return result;
}


/// Create a vector containing the imaginar part of the specified vector
/** 
 * 
 * 
 * @param x existing vector with complex values
 * 
 * @return pointer to new vector instance containing imaginary part of @c x
 */
inline Vector *imaginary(const Vector& x)
{
  Vector *result(x.clone());
  result->imaginary();
  return result;
}

/// Create a vector containing the complex conjugate of @c x
/** 
 * 
 * 
 * @param x existing vector with complex values
 * 
 * @return pointer to new vector instance containing the complex conjugate of @c x
 */
inline Vector* conjugate(const Vector& x)
{
  Vector *result(x.clone());
  result->conjugate();
  return result;
}



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
inline void add(const Vector& A, const Vector& B, Vector& result)
{
  result.equate(A);
  result.add(B);
}



} // namespace utility
} // namespace gridpack

#endif
