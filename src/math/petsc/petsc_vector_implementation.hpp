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
 * @date   2015-07-24 10:32:00 d3g096
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
#include "petsc/petsc_types.hpp"

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

  /// A flag to denote whether the library can be used directly
  /**
   * Some operations can be passed directly to the underlying library
   * if the TheType is the same as the PETSc type @e or the PETSc type
   * is complex.  This type computes and stores that flag. 
   * 
   */
  static const bool useLibrary = UsePetscLibrary<TheType>::value;

  /// The number of library elements used to represent a single vector element
  static const unsigned int elementSize = PetscElementSize<TheType>::value;

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
    : VectorImplementation<T>(comm), p_vwrap(comm, local_length*elementSize)
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

  /// Apply a specific unary operation to the vector
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

  /// Apply a specificy accumulator operation to the vector
  RealType p_applyAccumulator(base_accumulator_function<TheType, RealType>& op) const
  {
    RealType result;
    PetscErrorCode ierr;
    const Vec *v = p_vwrap.getVector();
    const PetscScalar *p;
    PetscInt n;
    ierr = VecGetLocalSize(*v, &n); CHKERRXX(ierr);
    ierr = VecGetArrayRead(*v, &p);  CHKERRXX(ierr);
    accumulator_operation<TheType, PetscScalar>(static_cast<unsigned int>(n), p, op);
    ierr = VecRestoreArrayRead(*v, &p); CHKERRXX(ierr);
    result = op.result();
    return result;
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

  void p_setOrAddElement(const IdxType& i, const TheType& x, InsertMode mode)
  {
    PetscErrorCode ierr;
    if (useLibrary) {
      try {
        Vec *v = p_vwrap.getVector();
        PetscScalar pv = equate<PetscScalar, TheType>(x);
        PetscInt pidx(i);
        ierr = VecSetValue(*v, pidx, pv, mode); CHKERRXX(ierr);
      } catch (const PETSC_EXCEPTION_TYPE& e) {
        throw PETScException(ierr, e);
      }
    } else {
      p_setOrAddElements(1, &i, &x, mode);
    }
  }

  /// Set an several elements (specialized)
  void p_setOrAddElements(const IdxType& n, const IdxType *i, const TheType *x, 
                          InsertMode mode)
  {
    PetscErrorCode ierr;
    try {
      Vec *v = p_vwrap.getVector();
      const IdxType *theindexes;
      PetscInt *idx = NULL;

      // try to avoid allocating an array for indexes: assumes that PetscInt == IdxType
      if (elementSize == 1) {
        theindexes = i;         
      } else {
        idx = new PetscInt[n*elementSize];
        for (int j = 0; j < n; ++j) {
          idx[j*elementSize] = i[j]*elementSize;
          if (elementSize > 1) 
            idx[j*elementSize+1] = idx[j*elementSize]+1;
        }
        theindexes = idx;
      }
      ValueTransferToLibrary<TheType, PetscScalar> trans(n, const_cast<TheType *>(x));
      trans.go();
      ierr = VecSetValues(*v, n*elementSize, &theindexes[0], trans.to(), mode); CHKERRXX(ierr);
      if (idx != NULL) {
        delete [] idx;
      }
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }


  /// Set an individual element (specialized)
  void p_setElement(const IdxType& i, const TheType& x)
  {
    p_setOrAddElement(i, x, INSERT_VALUES);
  }

  /// Set an several elements (specialized)
  void p_setElements(const IdxType& n, const IdxType *i, const TheType *x)
  {
    p_setOrAddElements(n, i, x, INSERT_VALUES);
  }

  /// Add to an individual element (specialized)
  void p_addElement(const IdxType& i, const TheType& x)
  {
    p_setOrAddElement(i, x, ADD_VALUES);
  }

  /// Add to an several elements (specialized)
  void p_addElements(const IdxType& n, const IdxType *i, const TheType *x)
  {
    p_setOrAddElements(n, i, x, ADD_VALUES);
  }

  /// Get an individual (local) element (specialized)
  void p_getElement(const IdxType& i, TheType& x) const
  {
    PetscErrorCode ierr;
    if (useLibrary) {
      try {
        const Vec *v = p_vwrap.getVector();
        PetscScalar pv;
        PetscInt pidx(i);
        ierr = VecGetValues(*v, 1, &pidx, &pv); CHKERRXX(ierr);
        x = equate<TheType, PetscScalar>(pv);
      } catch (const PETSC_EXCEPTION_TYPE& e) {
        throw PETScException(ierr, e);
      }
    } else {
      p_getElements(1, &i, &x);
    }
  }

  /// Get an several (local) elements (specialized)
  void p_getElements(const IdxType& n, const IdxType *i, TheType *x) const
  {
    PetscErrorCode ierr;
    try {
      const Vec *v = p_vwrap.getVector();
      std::vector<PetscScalar> px(n*elementSize);
      std::vector<PetscInt> idx(n*elementSize);
      for (int j = 0; j < n; ++j) {
        idx[j*elementSize] = i[j]*elementSize;
        if (elementSize > 1) 
          idx[j*elementSize+1] = idx[j*elementSize]+1;
      }
      ierr = VecGetValues(*v, n*elementSize, &idx[0], &px[0]); CHKERRXX(ierr);
      ValueTransferFromLibrary<PetscScalar, TheType> trans(n*elementSize, &px[0], x);
      trans.go();
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Get all of vector elements (on all processes)
  void p_getAllElements(TheType *x) const
  {
    unsigned int n(p_vwrap.size());
    std::vector<PetscScalar> px(n);
    p_vwrap.getAllElements(&px[0]);
    ValueTransferFromLibrary<PetscScalar, TheType> trans(n, &px[0], x);
    trans.go();
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
      if (useLibrary) {
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
    if (useLibrary) {
      Vec *vec = p_vwrap.getVector();
      PetscErrorCode ierr(0);
      try {
        PetscScalar px =
          gridpack::math::equate<PetscScalar, TheType>(x);
        ierr = VecScale(*vec, px); CHKERRXX(ierr);
      } catch (const PETSC_EXCEPTION_TYPE& e) {
        throw PETScException(ierr, e);
      }
    } else {
      gridpack::math::multiplyvalue<TheType> op(x);
      p_applyOperation(op);
    } 
  }

  /// Compute the vector L1 norm (sum of absolute value) (specialized)
  double p_norm1(void) const
  {
    double result;
    if (useLibrary) {
      result = p_vwrap.norm1();
    } else {
      l1_norm<TheType> op;
      double lresult(p_applyAccumulator(op));
      boost::mpi::all_reduce(this->communicator(), lresult, result, std::plus<double>());
    }
    return result;
  }

  /// Compute the vector L2 norm (root of sum of squares) (specialized)
  double p_norm2(void) const
  {
    double result;
    if (useLibrary) {
      result = p_vwrap.norm2();
    } else {
      l2_norm<TheType> op;
      double lresult(p_applyAccumulator(op));
      boost::mpi::all_reduce(this->communicator(), lresult, result, std::plus<double>());
      result = sqrt(result);
    }
    return result;
  }

  /// Compute the vector infinity (or maximum) norm (specialized)
  double p_normInfinity(void) const
  {
    double result;
    if (useLibrary) {
      result = p_vwrap.normInfinity();
    } else {
      infinity_norm<TheType> op;
      double lresult(p_applyAccumulator(op));
      boost::mpi::all_reduce(this->communicator(), lresult, result, boost::mpi::maximum<double>());
    }
    return result;
  }

  /// Replace all elements with its absolute value (specialized) 
  void p_abs(void)
  {
    if (useLibrary) {
      p_vwrap.abs();
    } else {
      absvalue<TheType> op;
      p_applyOperation(op);
    }
  }

  /// Replace all elements with their complex conjugate
  void p_conjugate(void)
  {
    if (useLibrary) {
      p_vwrap.conjugate();
    } else {
      conjugate_value<TheType> op;
      p_applyOperation(op);
    }
  }

  /// Replace all elements with its exponential (specialized)
  void p_exp(void)
  {
    if (useLibrary) {
      p_vwrap.exp();
    } else {
      exponential<TheType> op;
      p_applyOperation(op);
    }
  }

  /// Replace all elements with its reciprocal (specialized)
  void p_reciprocal(void)
  {
    if (useLibrary) {
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
