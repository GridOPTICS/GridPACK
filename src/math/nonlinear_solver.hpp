// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver.hpp
 * @author William A. Perkins
 * @date   2013-09-06 11:16:55 d3g096
 * 
 * @brief  
 * 
 * 
 */

#ifndef _nonlinear_solver_hpp_
#define _nonlinear_solver_hpp_

#include <boost/scoped_ptr.hpp>
#include <gridpack/math/nonlinear_solver_implementation.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolver
// -------------------------------------------------------------
class NonlinearSolver 
  : public parallel::WrappedDistributed,
    private utility::Uncopyable 
{
public:

  /// Default constructor.
  NonlinearSolver(const parallel::Communicator& comm,
                  const int& local_size,
                  JacobianBuilder form_jacobian,
                  FunctionBuilder form_function);

  /// Destructor
  ~NonlinearSolver(void);

  /// Solve w/ the specified initial estimated, put result in same vector
  void solve(Vector& x)
  {
    p_impl->solve(x);
  }

protected:

  /// Where things really happen
  boost::scoped_ptr<NonlinearSolverImplementation> p_impl;
  
};


} // namespace math
} // namespace gridpack


#endif
