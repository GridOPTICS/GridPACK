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
//  class PETScSparseParallelMatrixImplementation
// -------------------------------------------------------------
class PETScSparseParallelMatrixImplementation 
  : public MatrixImplementation
{
public:

  /// Default constructor.
  PETScSparseParallelMatrixImplementation(const parallel::Distribution& dist,
                                          const int& rows, const int& cols);

  /// Destructor
  ~PETScSparseParallelMatrixImplementation(void);

protected:

  /// The PETSc matrix representation
  Mat matrix_;
};

// -------------------------------------------------------------
//  class PETScSparseSerialMatrixImplementation
// -------------------------------------------------------------
class PETScSparseSerialMatrixImplementation 
  : public MatrixImplementation
{
public:

  /// Default constructor.
  PETScSparseSerialMatrixImplementation(const parallel::Distribution& dist,
                                        const int& rows, const int& cols);

  /// Destructor
  ~PETScSparseSerialMatrixImplementation(void);

protected:
  
  /// The PETSc matrix representation
  Mat matrix_;

};

} // namespace utility
} // namespace gridpack



#endif
