// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   linear_powerflow_solver.hpp
 * @author William A. Perkins
 * @date   2013-04-29 07:25:43 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _linear_powerflow_solver_hpp_
#define _linear_powerflow_solver_hpp_



#include "gridpack/applications/powerflow/powerflow_solver.hpp"
#include "gridpack/math/linear_solver.hpp"

namespace tracker {

// -------------------------------------------------------------
//  class LinearPowerflowSolver
// -------------------------------------------------------------
/// Solve power flow equations with a general linear system solver
/**
 * This is simply an interface to gridpack::math::LinearSolver.  The
 * power flow equations are simply handed to the linear solver.
 * Subclasses may manipulate the power flow equations before handing
 * them to the solver.
 * 
 */
class LinearPowerflowSolver : public PowerflowSolver 
{
public:

  /// Default constructor.
  LinearPowerflowSolver(math::Matrix *admittance, 
                        math::Vector *power_injection, 
                        math::Vector *bus_voltage);

  /// Destructor
  ~LinearPowerflowSolver(void);

protected:

  boost::scoped_ptr<math::LinearSolver> linear_solver_;

  /// Performs initial set system setup
  void solve_setup_(void) 
  {
    linear_solver_.reset(new math::LinearSolver(*admittance_));
    if (!linear_solver_->sane()) {
      throw gridpack::runtime_error("bad linear solver");
    }
  }

  /// Perform the solution
  void solve_(void)
  {
    linear_solver_->solve(*power_injection_, *bus_voltage_);
  }
};



} // namespace tracker


#endif
