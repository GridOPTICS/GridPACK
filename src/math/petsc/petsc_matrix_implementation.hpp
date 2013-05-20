// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_implementation.h
 * @author William A. Perkins
 * @date   2013-05-17 13:35:07 d3g096
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
#include "gridpack/math/matrix_implementation.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScMatrixImplementation
// -------------------------------------------------------------
class PETScMatrixImplementation
  : public MatrixImplementation
{
public:

  /// Default constructor
  PETScMatrixImplementation(const parallel::Communicator& comm,
                            const int& local_rows, const int& cols,
                            const bool& dense = false);

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
  void p_local_row_range(int& lo, int& hi) const;

  /// Get the total number of rows in this matrix (specialized)
  int p_rows(void) const;

  /// Get the number of local rows in this matirx (specialized)
  int p_local_rows(void) const;

  /// Get the number of columns in this matrix (specialized)
  int p_cols(void) const;

  // /// Set an individual element
  // void p_set_element(const int& i, const int& j, const complex_type& x);

  // /// Set an several element
  // void p_set_elements(const int *i, const int *j, const complex_type *x);

  // /// Set all elements in a row
  // void p_set_row(const int& i, const int *j, const complex_type *x);

  // /// Set all elements in a region
  // void p_set_region(const int& ni, const int& nj, 
  //                          const int *i, const int *j, const complex_type *x);

  // /// Add to  an individual element
  // void p_add_element(const int& i, const int& j, const complex_type& x);

  // /// Add to  an several element
  // void p_add_elements(const int *i, const int *j, const complex_type *x);

  // /// Add to  all elements in a row
  // void p_add_row(const int& i, const int *j, const complex_type *x);

  // /// Get an individual element
  // void p_get_element(const int& i, const int& j, complex_type& x) const;

  // /// Get an several element
  // void p_get_elements(const int *i, const int *j, complex_type *x) const;

  // /// Get all elements in a row
  // void p_get_row(const int& i, const int *j, complex_type *x) const;

  // /// Get all elements in a region
  // void p_get_region(const int& ni, const int& nj, 
  //                          const int *i, const int *j, complex_type *x) const;

  /// Make this instance ready to use
  void p_ready(void);

  /// Allow visits by implementation visitors
  void p_accept(ImplementationVisitor& visitor);

  /// Allow visits by implementation visitors
  void p_accept(ConstImplementationVisitor& visitor) const;

};

} // namespace utility
} // namespace gridpack



#endif
