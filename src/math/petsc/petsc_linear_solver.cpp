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
 * @date   2015-03-05 13:17:17 d3g096
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
template <typename T, typename I>
LinearSolverT<T, I>::LinearSolverT(LinearSolverT<T, I>::MatrixType& A)
  : parallel::WrappedDistributed(),
    utility::WrappedConfigurable(),
    utility::Uncopyable(),
    p_solver(new PETScLinearSolverImplementation<T, I>(A))
{
  p_setDistributed(p_solver.get());
  p_setConfigurable(p_solver.get());
  // empty
}

template
LinearSolverT<ComplexType>::LinearSolverT(LinearSolverT<ComplexType>::MatrixType& A);

template
LinearSolverT<RealType>::LinearSolverT(LinearSolverT<RealType>::MatrixType& A);

} // namespace math
} // namespace gridpack
