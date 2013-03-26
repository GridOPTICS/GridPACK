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
  PETScSparseParallelMatrixImplementation();

  /// Destructor
  ~PETScSparseParallelMatrixImplementation(void);

protected:
};

// -------------------------------------------------------------
//  class PETScSparseSerialMatrixImplementation
// -------------------------------------------------------------
class PETScSparseSerialMatrixImplementation 
  : public MatrixImplementation
{
public:

  /// Default constructor.
  PETScSparseSerialMatrixImplementation();

  /// Destructor
  ~PETScSparseSerialMatrixImplementation(void);

protected:
  
  

};

} // namespace utility
} // namespace gridpack



#endif
