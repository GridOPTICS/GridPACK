// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver.hpp
 * @author William A. Perkins
 * @date   2013-09-09 10:17:26 d3g096
 * 
 * @brief  
 * 
 * 
 */

#ifndef _nonlinear_solver_hpp_
#define _nonlinear_solver_hpp_

#include <boost/scoped_ptr.hpp>
#include "gridpack/math/nonlinear_solver_interface.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolver
// -------------------------------------------------------------
class NonlinearSolver 
  : public NonlinearSolverInterface
{
public:

  /// Default constructor.
  NonlinearSolver(const parallel::Communicator& comm,
                  const int& local_size,
                  JacobianBuilder form_jacobian,
                  FunctionBuilder form_function);

  /// Destructor
  ~NonlinearSolver(void);

};


} // namespace math
} // namespace gridpack


#endif
