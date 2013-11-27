// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver.cpp
 * @author William A. Perkins
 * @date   2013-11-08 11:49:58 d3g096
 * 
 * @brief  PETSc-specific implementation of NonlinearSolver
 * 
 * 
 */
// -------------------------------------------------------------

#include "nonlinear_solver.hpp"
#include "petsc/petsc_nonlinear_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolver
// -------------------------------------------------------------

// -------------------------------------------------------------
// NonlinearSolver:: constructors / destructor
// -------------------------------------------------------------
NonlinearSolver::NonlinearSolver(const parallel::Communicator& comm,
                                 const int& local_size,
                                 JacobianBuilder form_jacobian,
                                 FunctionBuilder form_function)
  : NonlinearSolverInterface()
{
  p_setImpl(new PetscNonlinearSolverImplementation(comm, local_size, 
                                                    form_jacobian, form_function));
  p_setDistributed(p_impl.get());
}

} // namespace math
} // namespace gridpack
