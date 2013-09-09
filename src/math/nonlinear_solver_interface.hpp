// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_interface.hpp
 * @author William A. Perkins
 * @date   2013-09-09 10:13:13 d3g096
 * 
 * @brief  
 * 
 * 
 */

#ifndef _nonlinear_solver_interface_hpp_
#define _nonlinear_solver_interface_hpp_

#include <boost/scoped_ptr.hpp>
#include <gridpack/math/nonlinear_solver_implementation.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolverInterface
// -------------------------------------------------------------
class NonlinearSolverInterface 
  : public parallel::WrappedDistributed,
    private utility::Uncopyable 
{
public:

  /// Default constructor.
  NonlinearSolverInterface();

  /// Destructor
  ~NonlinearSolverInterface(void);

  /// Solve w/ the specified initial estimated, put result in same vector
  void solve(Vector& x)
  {
    p_impl->solve(x);
  }

protected:

  /// Where things really happen
  boost::scoped_ptr<NonlinearSolverImplementation> p_impl;
  
  /// Set the implementation
  void p_set_impl(NonlinearSolverImplementation *impl);
};


} // namespace math
} // namespace gridpack


#endif
