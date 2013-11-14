// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

// -------------------------------------------------------------
/**
 * @file   petsc_dae_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2013-11-14 11:41:16 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_dae_solver_implementation_hpp_
#define _petsc_dae_solver_implementation_hpp_

#include <petscts.h>
#include "dae_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScDAESolverImplementation
// -------------------------------------------------------------
class PETScDAESolverImplementation 
  : public DAESolverImplementation
{
public:

  /// Default constructor.
  PETScDAESolverImplementation(const parallel::Communicator& comm, 
                               const int local_size,
                               DAEJacobianBuilder& jbuilder,
                               DAEFunctionBuilder& fbuilder);

  /// Destructor
  ~PETScDAESolverImplementation(void);

protected:

  /// The actual solver
  TS p_ts;

  /// The Jacobian matrix
  Mat *p_petsc_J;

  /// Do what is necessary to build this instance
  void p_build(const std::string& option_prefix);

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props);

  /// Initialize the system (specialized)
  void p_initialize(const double& t0,
                    const double& deltat0,
                    Vector& x0);

  /// Solve the system
  void p_solve(double& maxtime,
               int& maxsteps);

  /// Routine to assemble Jacobian that is sent to PETSc
  static PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec x, Vec xdot, 
                                      PetscReal a, Mat *jac, Mat *B, 
                                      MatStructure *flag, void *dummy);

  /// Routine to assemble RHS that is sent to PETSc
  static PetscErrorCode FormIFunction(TS ts, PetscReal t, Vec x, Vec xdot, 
                                      Vec F, void *dummy);
};



} // namespace math
} // namespace gridpack

#endif
