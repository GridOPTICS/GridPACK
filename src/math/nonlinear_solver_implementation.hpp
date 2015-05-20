// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2015-03-26 14:27:58 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------

#ifndef _nonlinear_solver_implementation_hpp_
#define _nonlinear_solver_implementation_hpp_

#include <boost/shared_ptr.hpp>
#include <gridpack/utilities/complex.hpp>
#include <gridpack/math/nonlinear_solver_functions.hpp>
#include <gridpack/math/nonlinear_solver_interface.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/configuration/configurable.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolverImplementation
// -------------------------------------------------------------
template <typename T, typename I = int>
class NonlinearSolverImplementation 
  : public NonlinearSolverInterface<T, I>, 
    public parallel::Distributed,
    public utility::Configurable,
    private utility::Uncopyable
{
public:

  typedef typename NonlinearSolverInterface<T, I>::VectorType VectorType;
  typedef typename NonlinearSolverInterface<T, I>::MatrixType MatrixType;
  typedef typename NLSBuilder<T, I>::Jacobian JacobianBuilder;
  typedef typename NLSBuilder<T, I>::Function FunctionBuilder;

  /// A functor to keep smart pointers from deleting their pointer
  struct null_deleter
  {
    void operator()(void const *) const { }
  };

  /// Default constructor.
  NonlinearSolverImplementation(const parallel::Communicator& comm,
                                const int& local_size,
                                JacobianBuilder form_jacobian,
                                FunctionBuilder form_function)
    : NonlinearSolverInterface<T, I>(),
      parallel::Distributed(comm), 
      utility::Configurable("NonlinearSolver"), 
      utility::Uncopyable(),
      p_J(), p_F(), 
      p_X((VectorType *)NULL, null_deleter()),  // pointer set by solve()
      p_jacobian(form_jacobian), 
      p_function(form_function),
      p_solutionTolerance(1.0e-05),
    p_functionTolerance(1.0e-10),
    p_maxIterations(50)
  {
    p_F.reset(new VectorType(this->communicator(), local_size));
    // std::cout << this->processor_rank() << ": "
    //           << "NonlinearSolverImplementation: construct Jacobian matrix: "
    //           << local_size << " x " << cols
    //           << std::endl;
    p_J.reset(new MatrixType(this->communicator(), local_size, local_size, Sparse));
  }

  /// Constructor that uses an existing Jacobian Matrix.
  NonlinearSolverImplementation(MatrixType& J,
                                JacobianBuilder form_jacobian,
                                FunctionBuilder form_function)
    : NonlinearSolverInterface<T, I>(),
      parallel::Distributed(J.communicator()), 
      utility::Configurable("NonlinearSolver"), 
      utility::Uncopyable(),
      p_J(&J, null_deleter()), p_F(), 
      p_X((VectorType *)NULL, null_deleter()),  // pointer set by solve()
      p_jacobian(form_jacobian), 
      p_function(form_function),
      p_solutionTolerance(1.0e-05),
      p_functionTolerance(1.0e-10),
      p_maxIterations(50)
  {
    p_F.reset(new VectorType(this->communicator(), J.localRows()));
  }


  /// Destructor
  virtual ~NonlinearSolverImplementation(void)
  {
    // empty
  }


protected:

  /// A place to compute and store the Jacobian
  boost::shared_ptr<MatrixType> p_J;

  /// A place to compute and store the RHS
  boost::shared_ptr<VectorType> p_F;

  /// A vector containing to current solution estimate
  boost::shared_ptr<VectorType> p_X;

  /// A thing to build a Jacobian
  JacobianBuilder p_jacobian;

  /// A thing to build the RHS
  FunctionBuilder p_function;

  /// The solution residual norm tolerance
  double p_solutionTolerance;

  /// The RHS function norm tolerance
  double p_functionTolerance;

  /// The maximum number of iterations to perform
  int p_maxIterations;

  /// Get the solution tolerance (specialized)
  double p_tolerance(void) const
  {
    return p_solutionTolerance;
  }

  /// Set the solver tolerance (specialized)
  void p_tolerance(const double& tol) 
  {
    p_solutionTolerance = tol;
  }

  /// Get the maximum iterations (specialized)
  int p_maximumIterations(void) const
  {
    return p_maxIterations;
  }

  /// Set the maximum iterations (specialized)
  void p_maximumIterations(const int& n) 
  {
    p_maxIterations = n;
  }

  /// Solve w/ using the specified initial guess, put solution in same vector
  void p_solve(VectorType& x)
  {
    // children should call this is their own p_solve()
    p_X.reset(&x, null_deleter());
  }

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props)
  {
    if (props) {
      p_solutionTolerance = props->get("SolutionTolerance", p_solutionTolerance);
      p_functionTolerance = props->get("FunctionTolerance", p_functionTolerance);
      p_maxIterations = props->get("MaxIterations", p_maxIterations);
    }
  }

};


} // namespace math
} // namespace gridpack


#endif
