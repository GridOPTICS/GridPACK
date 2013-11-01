// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   newton_raphson_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-10-29 07:47:51 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include "newton_raphson_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NewtonRaphsonSolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// NewtonRaphsonSolverImplementation:: constructors / destructor
// -------------------------------------------------------------
NewtonRaphsonSolverImplementation::NewtonRaphsonSolverImplementation(const parallel::Communicator& comm,
                                                                     const int& local_size,
                                                                     JacobianBuilder form_jacobian,
                                                                     FunctionBuilder form_function)
  : NonlinearSolverImplementation(comm, local_size, form_jacobian, form_function),
    p_tolerance(1.0e-03),
    p_max_iterations(1),
    p_linear_solver()
{
  this->configurationKey("NewtonRaphsonSolver");
}

NewtonRaphsonSolverImplementation::~NewtonRaphsonSolverImplementation(void)
{
}

// -------------------------------------------------------------
// NewtonRaphsonSolverImplementation::p_solve
// -------------------------------------------------------------
void
NewtonRaphsonSolverImplementation::p_solve(void)
{
  ComplexType stol(1.0e+30);
  ComplexType ftol(1.0e+30);
  int iter(0);

  boost::scoped_ptr<Vector> deltaX(p_X->clone());
  while (real(stol) > p_tolerance && iter < p_max_iterations) {
    p_function(*p_X, *p_F);
    p_F->scale(-1.0);
    p_jacobian(*p_X, *p_J);
    if (!p_linear_solver) {
      p_linear_solver.reset(new LinearSolver(*p_J));
      p_linear_solver->configure(this->p_configCursor);
    } else {
      p_linear_solver->setMatrix(*p_J);
    }
    deltaX->zero();
    p_linear_solver->solve(*p_F, *deltaX);
    stol = deltaX->norm2();
    ftol = p_F->norm2();
    p_X->add(*deltaX);
    iter += 1;
    if (this->processor_rank() == 0) {
      std::cout << "Newton-Raphson "
                << "iteration " << iter << ": "
                << "solution residual norm = " << real(stol) << ", "
                << "function norm = " << real(ftol)
                << std::endl;
    }
  }
}

// -------------------------------------------------------------
// NewtonRaphsonSolverImplementation::p_configure
// -------------------------------------------------------------
void
NewtonRaphsonSolverImplementation::p_configure(utility::Configuration::Cursor *props)
{
  if (props != NULL) {
    p_tolerance = props->get("Tolerance", p_tolerance);
    p_max_iterations = props->get("MaxIterations", p_max_iterations);
  }
}


} // namespace math
} // namespace gridpack
