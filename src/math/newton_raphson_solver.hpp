// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   newton_raphson_solver.hpp
 * @author William A. Perkins
 * @date   2013-09-09 10:24:25 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _newton_raphson_solver_hpp_
#define _newton_raphson_solver_hpp_

#include "gridpack/math/nonlinear_solver_interface.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NewtonRaphsonSolver
// -------------------------------------------------------------
class NewtonRaphsonSolver 
  : public NonlinearSolverInterface {
public:

  /// Default constructor.
  NewtonRaphsonSolver(const parallel::Communicator& comm,
                      const int& local_size,
                      JacobianBuilder form_jacobian,
                      FunctionBuilder form_function);
  
  /// Destructor
  ~NewtonRaphsonSolver(void);
};


} // namespace math
} // namespace gridpack

#endif
