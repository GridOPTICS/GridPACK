// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   petsc_linear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2013-10-08 10:08:46 d3g096
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
  PETScLinearSolverImplementation(const Matrix& A);

  /// Destructor
  ~PETScLinearSolverImplementation(void);

protected:

  /// The PETSc linear solver 
  KSP p_KSP;

  /// Do what is necessary to build this instance
  void p_build(const std::string& option_prefix);

  /// Solve w/ the specified RHS and estimate (result in x)
  void p_solve(const Vector& b, Vector& x) const;

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::Cursor *props);

  /// Use different coefficient matrix (or A w/ new values) (specialized)
  void p_setMatrix(void);

  /// Allow visits by implementation visitors
  void p_accept(ImplementationVisitor& visitor);

  /// Allow visits by implementation visitors
  void p_accept(ConstImplementationVisitor& visitor) const;

};

} // namespace math
} // namespace gridpack

#endif
