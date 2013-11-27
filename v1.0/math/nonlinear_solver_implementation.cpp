// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-11-08 09:04:45 d3g096
 * 
 * @brief  Abstract class NonlinearSolverImplementation implementation 
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/mpi/collectives.hpp>
#include "nonlinear_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
// struct null_deleter
// -------------------------------------------------------------
struct null_deleter
{
  void operator()(void const *) const { }
};

// -------------------------------------------------------------
//  class NonlinearSolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// NonlinearSolverImplementation:: constructors / destructor
// -------------------------------------------------------------
NonlinearSolverImplementation::NonlinearSolverImplementation(const parallel::Communicator& comm,
                                                             const int& local_size,
                                                             JacobianBuilder form_jacobian,
                                                             FunctionBuilder form_function)
  : parallel::Distributed(comm), 
    utility::Configurable("NonlinearSolver"), 
    utility::Uncopyable(),
    p_J(), p_F(), 
    p_X((Vector *)NULL, null_deleter()),  // pointer set by solve()
    p_jacobian(form_jacobian), 
    p_function(form_function),
    p_solutionTolerance(1.0e-05),
    p_functionTolerance(1.0e-10),
    p_maxIterations(50)
{
  p_F.reset(new Vector(this->communicator(), local_size));
  int cols;
  all_reduce(this->communicator(), local_size, cols, std::plus<int>());
  // std::cout << this->processor_rank() << ": "
  //           << "NonlinearSolverImplementation: construct Jacobian matrix: "
  //           << local_size << " x " << cols
  //           << std::endl;
  p_J.reset(new Matrix(this->communicator(), local_size, cols, Matrix::Sparse));
}

NonlinearSolverImplementation::~NonlinearSolverImplementation(void)
{
  // empty
}

// -------------------------------------------------------------
// NonlinearSolverImplementation::solve
// -------------------------------------------------------------
void 
NonlinearSolverImplementation::solve(Vector &x) 
{
  p_X.reset(&x, null_deleter());
  this->p_solve();
}

// -------------------------------------------------------------
// NonlinearSolverImplementation::p_configure
// -------------------------------------------------------------
void
NonlinearSolverImplementation::p_configure(utility::Configuration::CursorPtr props)
{
  if (props) {
    p_solutionTolerance = props->get("SolutionTolerance", p_solutionTolerance);
    p_functionTolerance = props->get("FunctionTolerance", p_functionTolerance);
    p_maxIterations = props->get("MaxIterations", p_maxIterations);
  }
}


} // namespace math
} // namespace gridpack
