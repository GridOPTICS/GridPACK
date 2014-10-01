// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   linear_matrix_solver.cpp
 * @author William A. Perkins
 * @date   2013-11-08 11:49:59 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "linear_matrix_solver.hpp"
#include "basic_linear_matrix_solver_implementation.hpp"
#include "petsc_linear_matrix_solver_impl.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearMatrixSolver
// -------------------------------------------------------------

// -------------------------------------------------------------
// LinearMatrixSolver:: constructors / destructor
// -------------------------------------------------------------
LinearMatrixSolver::LinearMatrixSolver(const Matrix& A)
  : parallel::WrappedDistributed(),
    utility::WrappedConfigurable(),
    utility::Uncopyable()
{
#ifdef PETSC_HAVE_SUPERLU
  p_impl.reset(new PetscLinearMatrixSolverImplementation(A));
#else
  p_impl.reset(new BasicLinearMatrixSolverImplementation(A));
#endif
  p_setDistributed(p_impl.get());
  p_setConfigurable(p_impl.get());
}


} // namespace math
} // namespace gridpack
