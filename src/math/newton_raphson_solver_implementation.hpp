// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: newton_raphson_solver_implementation.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created September  9, 2013 by William A. Perkins
// Last Change: 2013-09-09 13:01:42 d3g096
// -------------------------------------------------------------


#ifndef _newton_raphson_solver_implementation_hpp_
#define _newton_raphson_solver_implementation_hpp_

#include <boost/scoped_ptr.hpp>
#include "nonlinear_solver_implementation.hpp"
#include "nonlinear_solver_functions.hpp"
#include "linear_solver.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NewtonRaphsonSolverImplementation
// -------------------------------------------------------------
class NewtonRaphsonSolverImplementation 
  : public NonlinearSolverImplementation
{
public:

  /// Default constructor.
  NewtonRaphsonSolverImplementation(const parallel::Communicator& comm,
                                    const int& local_size,
                                    JacobianBuilder form_jacobian,
                                    FunctionBuilder form_function);

  /// Destructor
  ~NewtonRaphsonSolverImplementation(void);

  /// Get the solver tolerance
  double tolerance(void) const
  {
    return p_tolerance;
  }

  /// Set the solver tolerance
  void tolerance(const double& tol) 
  {
    p_tolerance = tol;
  }

  /// Get the maximum iterations
  int maximum_iterations(void) const
  {
    return p_max_iterations;
  }

  /// Set the maximum iterations
  void maximum_iterations(const int& n) 
  {
    p_max_iterations = n;
  }

protected:

  /// Solver tolerance goal
  double p_tolerance;

  /// Maximum number of iterations to perform
  int p_max_iterations;

  /// The linear solver
  boost::scoped_ptr<LinearSolver> p_linear_solver;

  /// Solve w/ using the specified initial guess (specialized)
  void p_solve(void);

};



} // namespace math
} // namespace gridpack
#endif
