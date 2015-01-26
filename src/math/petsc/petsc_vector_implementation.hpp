// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_vector_implementation.hpp
 * @author William A. Perkins
 * @date   2015-01-26 09:55:58 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_vector_implementation_h_
#define _petsc_vector_implementation_h_

#include "complex_operators.hpp"
#include "value_transfer.hpp"
#include "vector_implementation.hpp"
#include "petsc/petsc_vector_wrapper.hpp"
#include "petsc/petsc_exception.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScVectorImplementation
// -------------------------------------------------------------
/// Vector implementation based on the PETSc library
/**
 * 
 * 
 */
template <typename T, typename I = int>
class PETScVectorImplementation 
  : public VectorImplementation<T, I> {
public:

  typedef typename VectorImplementation<T, I>::IdxType IdxType;
  typedef typename VectorImplementation<T, I>::TheType TheType;

  /// A flag (type) to denote whether the library can be used
  /**
   * Some operations can be passed directly to the underlying library
   * if the TheType is the same as the PETSc type @e or the PETSc type
   * is complex.  This type computes and stores that flag. 
   * 
   */

  typedef 
  typename boost::mpl::bool_<
    boost::is_same<TheType, PetscScalar>::value ||
    boost::is_same<ComplexType, PetscScalar>::value
  >::type useLibrary;

  /// The number of library elements used to represent a single vector element
  static const unsigned int elementSize = 
    boost::mpl::if_< useLibrary, 
                     boost::mpl::int_<1>, 
                     boost::mpl::int_<2> >::type::value;

    // storage_size<TheType, PetscScalar>::value;



  /// Default constructor.
  /** 
   * @e Collective on @c comm.
   * 
   * @param comm parallel environment
   * @param local_length number of entries to be locally owned
   * 
   * @return 
   */
  PETScVectorImplementation(const parallel::Communicator& comm,
                            const IdxType& local_length)
    : VectorImplementation<T>(comm), p_vwrap(comm, local_length)
  { }

  /// Construct from an existing PETSc vector
  PETScVectorImplementation(Vec& pvec, const bool& copyvec = true)
    : VectorImplementation<T>(PetscVectorWrapper::getCommunicator(pvec)), 
      p_vwrap(pvec, copyvec)
  { }

  /// Destructor
  /** 
   * @e Collective
   * 
   */
  ~PETScVectorImplementation(void) {}

protected:

  /// A vector of TheType
  typedef std::vector<TheType> TheVector;

  /// Where the actual vector is stored
  PetscVectorWrapper p_vwrap;

  // -------------------------------------------------------------
  // p_applyOperation
  // -------------------------------------------------------------
  void p_applyOperation(base_unary_function<TheType>& op)
  {
    PetscErrorCode ierr;
    Vec *v = p_vwrap.getVector();
    PetscScalar *p;
    PetscInt n;
    ierr = VecGetLocalSize(*v, &n); CHKERRXX(ierr);
    ierr = VecGetArray(*v, &p);  CHKERRXX(ierr);
    unary_operation<TheType, PetscScalar>(static_cast<unsigned int>(n), 
                                          p, op);
    ierr = VecRestoreArray(*v, &p); CHKERRXX(ierr);
    this->ready();
  }

  /// Get the global vector length
  IdxType p_size(void) const
  {
    return p_vwrap.size()/elementSize;
  }

  /// Get the size of the vector local part
  IdxType p_localSize(void) const
  {
    return p_vwrap.localSize()/elementSize;
  }

  /// Get the local min/max global indexes (specialized)
  void p_localIndexRange(IdxType& lo, IdxType& hi) const
  {
    PetscInt plo, phi;
    p_vwrap.localIndexRange(plo, phi);
    lo = plo/elementSize;
    hi = phi/elementSize;
  }

  /// Set an individual element (specialized)
  /** 
   * If you try to set an off processor value, it will be ignored
   * 
   * @param i global vector index
   * @param x value
   */
  void p_setElement(const IdxType& i, const TheType& x)
  {
    PetscErrorCode ierr;
    try {
      Vec *v = p_vwrap.getVector();
      PetscScalar px(x);
      PetscInt idx(i/elementSize);
      ierr = VecSetValue(*v, idx, px, INSERT_VALUES); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Set an several elements (specialized)
  void p_setElements(const IdxType& n, const IdxType *i, const TheType *x)
  {
    // FIXME: should be able to set multiple values
    for (int idx = 0; idx < n; ++idx) {
      this->p_setElement(i[idx], x[idx]);
    }
  }

  /// Add to an individual element (specialized)
  void p_addElement(const IdxType& i, const TheType& x)
  {
    PetscErrorCode ierr;
    try {
      Vec *v = p_vwrap.getVector();
      PetscScalar px(x);
      ierr = VecSetValue(*v, i, px, ADD_VALUES); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Add to an several elements (specialized)
  void p_addElements(const IdxType& n, const IdxType *i, const TheType *x)
  {
    // FIXME: should be able to add multiple values
    for (int idx = 0; idx < n; ++idx) {
      this->p_addElement(i[idx], x[idx]);
    }
  }

  /// Get an individual (local) element (specialized)
  void p_getElement(const IdxType& i, TheType& x) const
  {
    PetscErrorCode ierr(0);
    try {
      const Vec *v = p_vwrap.getVector();
      PetscScalar y;
      PetscInt idx(i);
      ierr = VecGetValues(*v, 1, &idx, &y); CHKERRXX(ierr);
      x = equate<TheType, PetscScalar>(y);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Get an several (local) elements (specialized)
  void p_getElements(const IdxType& n, const IdxType *i, TheType *x) const
  {
    // FIXME: should be able to get multiple values
    for (int idx = 0; idx < n; ++idx) {
      this->p_getElement(i[idx], x[idx]);
    }
  }

  /// Get all of vector elements (on all processes)
  void p_getAllElements(TheType *x) const
  {
    std::vector<PetscScalar> px(this->size());
    p_vwrap.getAllElements(&px[0]);
    for (size_t i = 0; i < px.size(); ++i) {
      x[i] = equate<TheType, PetscScalar>(px[i]);
    }
  }

  /// Make all the elements zero (specialized)
  void p_zero(void)
  {
    p_vwrap.zero();
  }

  /// Make all the elements the specified value (specialized)
  void p_fill(const TheType& value)
  {
    PetscErrorCode ierr(0);
    try {
      if (useLibrary::value) {
        Vec *v = p_vwrap.getVector();
        PetscScalar pv = 
          equate<PetscScalar, TheType>(value);
        ierr = VecSet(*v, pv); CHKERRXX(ierr);
      } else {
        setvalue<TheType> op(value);
        p_applyOperation(op);
      }
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }  

  /// Scale all elements by a single value
  void p_scale(const TheType& x)
  {
    if (useLibrary::value) {
      Vec *vec = p_vwrap.getVector();
      PetscErrorCode ierr(0);
      try {
        PetscScalar px(x);
        ierr = VecScale(*vec, px); CHKERRXX(ierr);
      } catch (const PETSC_EXCEPTION_TYPE& e) {
        throw PETScException(ierr, e);
      }
    } else {
      gridpack::math::multiply<TheType> op(x);
      p_applyOperation(op);
    } 
  }

  /// Compute the vector L1 norm (sum of absolute value) (specialized)
  double p_norm1(void) const
  {
    double result;
    if (useLibrary::value) {
      result = p_vwrap.norm1();
    } else {
      BOOST_ASSERT(false);
    }
    return result;
  }

  /// Compute the vector L2 norm (root of sum of squares) (specialized)
  double p_norm2(void) const
  {
    double result;
    if (useLibrary::value) {
      result = p_vwrap.norm2();
    } else {
      BOOST_ASSERT(false);
    }
    return result;
  }

  /// Compute the vector infinity (or maximum) norm (specialized)
  double p_normInfinity(void) const
  {
    double result;
    if (useLibrary::value) {
      result = p_vwrap.normInfinity();
    } else {
      BOOST_ASSERT(false);
    }
    return result;
  }

  /// Replace all elements with its absolute value (specialized) 
  void p_abs(void)
  {
    if (useLibrary::value) {
      p_vwrap.abs();
    } else {
      absvalue<TheType> op;
      p_applyOperation(op);
    }
  }

  /// Replace all elements with their complex conjugate
  void p_conjugate(void)
  {
    if (useLibrary::value) {
      p_vwrap.conjugate();
    } else {
      BOOST_ASSERT(false);
    }
  }

  /// Replace all elements with its exponential (specialized)
  void p_exp(void)
  {
    if (useLibrary::value) {
      p_vwrap.exp();
    } else {
      exponential<TheType> op;
      p_applyOperation(op);
    }
  }

  /// Replace all elements with its reciprocal (specialized)
  void p_reciprocal(void)
  {
    if (useLibrary::value) {
      p_vwrap.reciprocal();
    } else {
      reciprocal<TheType> op;
      p_applyOperation(op);
    }
  }

  /// Make this instance ready to use
  void p_ready(void)
  {
    p_vwrap.ready();
  }

  /// Allow visits by implemetation visitor
  void p_accept(ImplementationVisitor& visitor)
  {
    p_vwrap.accept(visitor);
  }

  /// Print to named file or standard output (specialized)
  void p_print(const char* filename = NULL) const 
  { 
    p_vwrap.print(filename); 
  }

  /// Save, in MatLAB format, to named file (collective) (specialized)
  void p_save(const char *filename) const 
  { 
    p_vwrap.save(filename); 
  }

  /// Load from a named file of whatever binary format the math library uses (specialized)
  void p_loadBinary(const char *filename) 
  { 
    p_vwrap.loadBinary(filename); 
  }

  /// Save to named file in whatever binary format the math library uses (specialized)
  void p_saveBinary(const char *filename) const 
  { 
    p_vwrap.saveBinary(filename); 
  }

  /// Allow visits by implemetation visitor
  void p_accept(ConstImplementationVisitor& visitor) const
  {
    p_vwrap.accept(visitor);
  }

  /// Make an exact replica of this instance (specialized)
  VectorImplementation<T> *p_clone(void) const
  {
    parallel::Communicator comm(this->communicator());
    IdxType local_size(this->localSize());
    
    PETScVectorImplementation<T> *result = 
      new PETScVectorImplementation<T>(comm, local_size);
    PetscErrorCode ierr;
    
    Vec *to_vec(result->p_vwrap.getVector());

    try {
      const Vec *v = p_vwrap.getVector();
      ierr = VecCopy(*v, *to_vec); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
    return result;
  }
};


} // namespace math
} // namespace gridpack

#endif
