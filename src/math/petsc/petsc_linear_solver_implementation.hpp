// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_linear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2013-10-09 13:24:08 d3g096
 * 
 * @brief  
 * 
 * 
 */
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
