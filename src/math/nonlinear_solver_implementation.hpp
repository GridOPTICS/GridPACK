// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2013-08-12 11:14:43 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------

#ifndef _nonlinear_solver_implementation_hpp_
#define _nonlinear_solver_implementation_hpp_

#include <boost/shared_ptr.hpp>
#include <gridpack/math/nonlinear_solver_functions.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolverImplementation
// -------------------------------------------------------------
class NonlinearSolverImplementation 
  : public parallel::Distributed,
    private utility::Uncopyable
{
public:

  /// Default constructor.
  NonlinearSolverImplementation(const parallel::Communicator& comm,
                                const int& local_size,
                                JacobianBuilder form_jacobian,
                                FunctionBuilder form_function);

  /// Destructor
  virtual ~NonlinearSolverImplementation(void);

  /// Solve w/ using the specified initial guess, put solution in same vector
  void solve(Vector& x);

protected:

  /// A place to compute and store the Jacobian
  boost::shared_ptr<Matrix> p_J;

  /// A place to compute and store the RHS
  boost::shared_ptr<Vector> p_F;

  /// A vector containing to current solution estimate
  boost::shared_ptr<Vector> p_X;

  /// A thing to build a Jacobian
  JacobianBuilder p_jacobian;

  /// A thing to build the RHS
  FunctionBuilder p_function;

  /// Solve w/ using the specified initial guess (specialized)
  virtual void p_solve(void) = 0;
};


} // namespace math
} // namespace gridpack


#endif
