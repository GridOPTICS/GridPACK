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
 * @date   2017-12-05 10:16:00 d3g096
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
template <typename T, typename I>
LinearMatrixSolverT<T, I>::LinearMatrixSolverT(LinearMatrixSolverT<T, I>::MatrixType& A)
  : BaseLinearMatrixSolverInterface<T, I>(),
    parallel::WrappedDistributed(),
    utility::WrappedConfigurable(),
    utility::Uncopyable()
{
#if defined(PETSC_HAVE_SUPERLU_DIST) || defined(PETSC_HAVE_MUMPS)
  p_impl.reset(new PetscLinearMatrixSolverImplementation<T, I>(A));
#else
  p_impl.reset(new BasicLinearMatrixSolverImplementation<T, I>(A));
#endif
  p_setDistributed(p_impl.get());
  p_setConfigurable(p_impl.get());
}


template
LinearMatrixSolverT<ComplexType>::LinearMatrixSolverT(LinearMatrixSolverT<ComplexType>::MatrixType& A);

template
LinearMatrixSolverT<RealType>::LinearMatrixSolverT(LinearMatrixSolverT<RealType>::MatrixType& A);

} // namespace math
} // namespace gridpack
