// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_implementation.h
 * @author William A. Perkins
 * @date   Mon Mar 25 12:26:00 2013
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
#ifndef _petsc_matrix_implementation_h_
#define _petsc_matrix_implementation_h_

#include <petscmat.h>
#include "matrix_implementation.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScMatrixImplementation
// -------------------------------------------------------------
class PETScMatrixImplementation
  : public MatrixImplementation
{
public:

  /// Default constructor.
  PETScMatrixImplementation(const parallel::Distribution& dist,
                            const int& rows, const int& cols);

  /// Destructor
  ~PETScMatrixImplementation(void);

protected:

  /// The PETSc matrix representation
  Mat matrix_;

  /// Set an individual element
  void set_element_(const int& i, const int& j, const complex_type& x);

  /// Set an several element
  void set_elements_(const int *i, const int *j, const complex_type *x);

  /// Set all elements in a row
  void set_row_(const int& i, const int *j, const complex_type *x);

  /// Set all elements in a region
  void set_region_(const int& ni, const int& nj, 
                   const int *i, const int *j, const complex_type *x) = 0;

  /// Add to  an individual element
  void add_element_(const int& i, const int& j, const complex_type& x);

  /// Add to  an several element
  void add_elements_(const int *i, const int *j, const complex_type *x);

  /// Add to  all elements in a row
  void add_row_(const int& i, const int *j, const complex_type *x);

  /// Get an individual element
  void get_element_(const int& i, const int& j, const complex_type& x);

  /// Get an several element
  void get_elements_(const int *i, const int *j, const complex_type *x);

  /// Get all elements in a row
  void get_row_(const int& i, const int *j, const complex_type *x);

  /// Get all elements in a region
  void get_region_(const int& ni, const int& nj, 
                   const int *i, const int *j, const complex_type *x);

};



// -------------------------------------------------------------
//  class PETScSparseParallelMatrixImplementation
// -------------------------------------------------------------
class PETScSparseParallelMatrixImplementation 
  : public PETScMatrixImplementation
{
public:

  /// Default constructor.
  PETScSparseParallelMatrixImplementation(const parallel::Distribution& dist,
                                          const int& rows, const int& cols);

  /// Destructor
  ~PETScSparseParallelMatrixImplementation(void);
};

// -------------------------------------------------------------
//  class PETScSparseSerialMatrixImplementation
// -------------------------------------------------------------
class PETScSparseSerialMatrixImplementation 
  : public PETScMatrixImplementation
{
public:

  /// Default constructor.
  PETScSparseSerialMatrixImplementation(const parallel::Distribution& dist,
                                        const int& rows, const int& cols);

  /// Destructor
  ~PETScSparseSerialMatrixImplementation(void);

};

// -------------------------------------------------------------
//  class PETScDenseParallelMatrixImplementation
// -------------------------------------------------------------
class PETScDenseParallelMatrixImplementation 
  : public PETScMatrixImplementation
{
public:

  /// Default constructor.
  PETScDenseParallelMatrixImplementation(const parallel::Distribution& dist,
                                          const int& rows, const int& cols);

  /// Destructor
  ~PETScDenseParallelMatrixImplementation(void);
};

// -------------------------------------------------------------
//  class PETScDenseSerialMatrixImplementation
// -------------------------------------------------------------
class PETScDenseSerialMatrixImplementation 
  : public PETScMatrixImplementation
{
public:

  /// Default constructor.
  PETScDenseSerialMatrixImplementation(const parallel::Distribution& dist,
                                        const int& rows, const int& cols);

  /// Destructor
  ~PETScDenseSerialMatrixImplementation(void);

};

} // namespace utility
} // namespace gridpack



#endif
