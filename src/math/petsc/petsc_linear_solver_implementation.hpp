// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   petsc_linear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   Mon Apr  1 09:08:49 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  1, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _petsc_linear_solver_implementation_hpp_
#define _petsc_linear_solver_implementation_hpp_

#include <petscksp.h>
#include "linear_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScLinearSolverImplementation
// -------------------------------------------------------------
class PETScLinearSolverImplementation 
  : public LinearSolverImplementation {
public:

  /// Default constructor.
  PETScLinearSolverImplementation(const parallel::Distribution& dist,
                                  const Matrix& A);

  /// Destructor
  ~PETScLinearSolverImplementation(void);

protected:

  /// The coefficient matrix in PETSc form
  Mat *A;

  /// The PETSc linear solver 
  KSP ksp;

  /// Solve w/ the specified RHS and estimate (result in x)
  void solve_(const Vector& b, Vector& x) const;

  /// Allow visits by implemetation visitor
  void accept_(ImplementationVisitor& visitor);
  {
    visitor->visit(*this);
  }

};

} // namespace math
} // namespace gridpack

#endif
