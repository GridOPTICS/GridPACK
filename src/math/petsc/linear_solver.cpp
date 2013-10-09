// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   linear_solver.cpp
 * @author William A. Perkins
 * @date   2013-10-09 13:22:20 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "linear_solver.hpp"
#include "petsc/petsc_linear_solver_implementation.hpp"

namespace gridpack {
namespace math {


// -------------------------------------------------------------
//  class LinearSolver
// -------------------------------------------------------------

// -------------------------------------------------------------
// LinearSolver:: constructors / destructor
// -------------------------------------------------------------
LinearSolver::LinearSolver(const Matrix& A)
  : parallel::WrappedDistributed(),
    utility::Uncopyable(),
    p_solver(new PETScLinearSolverImplementation(A))
{
  p_set_distributed(p_solver.get());
  // empty
}


} // namespace math
} // namespace gridpack
