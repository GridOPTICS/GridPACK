// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_implementation.h
 * @author William A. Perkins
 * @date   2013-10-09 13:24:43 d3g096
 * 
 * @brief  
 * 
 * 
 */
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
/// Matrix implementation based on the PETSc library
/**
 * 
 * 
 */
class PETScMatrixImplementation
  : public MatrixImplementation
{
public:

  /// Default constructor
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const int& local_rows, const int& cols,
                            const bool& dense = false);

  /// Make a new instance from an existing PETSc matrix
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const Mat& m);

  /// Destructor
  ~PETScMatrixImplementation(void);

  /// Get the PETSc matrix
  Mat *get_matrix(void)
  {
    return &p_matrix;
  }

  /// Get the PETSc matrix (const version)
  const Mat *get_matrix(void) const
  {
    return &p_matrix;
  }

protected:

  /// The PETSc matrix representation
  Mat p_matrix;

  /// Get the global index range of the locally owned rows (specialized)
  void p_localRowRange(int& lo, int& hi) const;

  /// Get the total number of rows in this matrix (specialized)
  int p_rows(void) const;

  /// Get the number of local rows in this matirx (specialized)
  int p_localRows(void) const;

  /// Get the number of columns in this matrix (specialized)
  int p_cols(void) const;

  /// Set an individual element
  void p_setElement(const int& i, const int& j, const ComplexType& x);

  /// Set an several element
  void p_setElements(const int& n, const int *i, const int *j, const ComplexType *x);

  // /// Set all elements in a row
  // void p_set_row(const int& i, const int *j, const ComplexType *x);

  // /// Set all elements in a region
  // void p_set_region(const int& ni, const int& nj, 
  //                          const int *i, const int *j, const ComplexType *x);

  /// Add to  an individual element
  void p_addElement(const int& i, const int& j, const ComplexType& x);

  /// Add to  an several element
  void p_addElements(const int& n, const int *i, const int *j, const ComplexType *x);

  // /// Add to  all elements in a row
  // void p_add_row(const int& i, const int *j, const ComplexType *x);

  /// Get an individual element
  void p_getElement(const int& i, const int& j, ComplexType& x) const;

  /// Get an several element
  void p_getElements(const int& n, const int *i, const int *j, ComplexType *x) const;

  // /// Get all elements in a row
  // void p_get_row(const int& i, const int *j, ComplexType *x) const;

  // /// Get all elements in a region
  // void p_get_region(const int& ni, const int& nj, 
  //                          const int *i, const int *j, ComplexType *x) const;

  /// Make this instance ready to use
  void p_ready(void);

  /// Allow visits by implementation visitors
  void p_accept(ImplementationVisitor& visitor);

  /// Allow visits by implementation visitors
  void p_accept(ConstImplementationVisitor& visitor) const;

  /// Make an exact replica of this instance (specialized)
  MatrixImplementation *p_clone(void) const;

};

} // namespace utility
} // namespace gridpack



#endif
