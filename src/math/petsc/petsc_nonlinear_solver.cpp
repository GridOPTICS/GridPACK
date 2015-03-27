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
 * @date   2015-03-26 12:54:41 d3g096
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
template <typename T, typename I>
NonlinearSolverT<T, I>::NonlinearSolverT(const parallel::Communicator& comm,
                                         const int& local_size,
                                         NonlinearSolverT<T, I>::JacobianBuilder form_jacobian,
                                         NonlinearSolverT<T, I>::FunctionBuilder form_function)
  : NonlinearSolverInterface<T, I>(),
    parallel::WrappedDistributed(),
    utility::WrappedConfigurable(),
    utility::Uncopyable()
{
  p_setImpl(new PetscNonlinearSolverImplementation<T, I>(comm, local_size, 
                                                         form_jacobian, form_function));
}

template
NonlinearSolverT<ComplexType>::NonlinearSolverT(const parallel::Communicator& comm,
                                                const int& local_size,
                                                NonlinearSolverT<ComplexType>::JacobianBuilder form_jacobian,
                                                NonlinearSolverT<ComplexType>::FunctionBuilder form_function);

template
NonlinearSolverT<RealType>::NonlinearSolverT(const parallel::Communicator& comm,
                                             const int& local_size,
                                             NonlinearSolverT<RealType>::JacobianBuilder form_jacobian,
                                             NonlinearSolverT<RealType>::FunctionBuilder form_function);

template <typename T, typename I>
NonlinearSolverT<T, I>::NonlinearSolverT(NonlinearSolverT<T, I>::MatrixType& J,
                                         NonlinearSolverT<T, I>::JacobianBuilder form_jacobian,
                                         NonlinearSolverT<T, I>::FunctionBuilder form_function)
  : NonlinearSolverInterface<T, I>(),
    parallel::WrappedDistributed(),
    utility::WrappedConfigurable(),
    utility::Uncopyable()
{
  p_setImpl(new PetscNonlinearSolverImplementation<T, I>(J, form_jacobian, form_function));
}

template
NonlinearSolverT<ComplexType>::NonlinearSolverT(NonlinearSolverT<ComplexType>::MatrixType& J,
                                                NonlinearSolverT<ComplexType>::JacobianBuilder form_jacobian,
                                                NonlinearSolverT<ComplexType>::FunctionBuilder form_function);

template
NonlinearSolverT<RealType>::NonlinearSolverT(NonlinearSolverT<RealType>::MatrixType& J,
                                             NonlinearSolverT<RealType>::JacobianBuilder form_jacobian,
                                             NonlinearSolverT<RealType>::FunctionBuilder form_function);

} // namespace math
} // namespace gridpack
