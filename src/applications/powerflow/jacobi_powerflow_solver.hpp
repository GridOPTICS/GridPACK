// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   jacobi_powerflow_solver.hpp
 * @author William A. Perkins
 * @date   2013-04-26 15:45:34 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------

#ifndef _jacobi_powerflow_solver_hpp_
#define _jacobi_powerflow_solver_hpp_

#include "gridpack/applications/powerflow/powerflow_solver.hpp"

namespace tracker {

// -------------------------------------------------------------
//  class JacobiPowerflowSolver
// -------------------------------------------------------------
/**
 * This implements the Jacobi iterative scheme (sometimes called the
 * Gauss method) for power flow equations.  This is only an example
 * for talking purposes.  You would never do this because there are
 * better low-level implementations.  The idea is that a variety of
 * methods can be subclassed from PowerflowSolver.
 *
 * 
 * 
 */

class JacobiPowerflowSolver 
  : public PowerflowSolver
{
protected:

  math::Matrix *R_;
  math::Vector *Dinv, *xnew;

  /// Performs initial set system setup
  void solve_setup_(void) 
  {
    // check to make sure the vector sizes match :: admittance_
    R_ = math::clone(*admittance_);
    
    Dinv_ = math::diagonal(*R_);
    Dinv_->reciprocal();

    Enew_ = math::clone(*bus_voltage_);
    Enew_->zero();

    R_->mutiply_diagonal(*Enew); // zeros diagonal
  }

  /// Solves the power flow equations
  /** 
   * This performs a specific number of Jacobi iterations (sometimes
   * called Gauss iterations) on the power flow equations.  
   * 
   */
  void solve_(void)
  {
    // these should probably be attributes filled @ construction
    int maxiter(this->get_param<int>("iterations"));
    int relax(this->get_param<math::complex_type>("relaxation");

    for (int i = 0; i < maxiter; ++i) {
      math::multiply(*R_, *bus_voltage_, *Enew);
      Enew_->subtract(power_injection_);
      Enew_->scale(-1.0);
      Enew_->element_multiply(Dinv_);
      Enew_->scale(relax);      // under/over relaxation
      bus_voltage_->copy(Enew_);
    }
  }

public:

  /// Default constructor.
  JacobiPowerflowSolver(math::Matrix *Ybus)
    : PowerflowSolver(Ybus)

  /// Destructor
  ~JacobiPowerflowSolver(void);
};


} // namespace tracker


#endif
