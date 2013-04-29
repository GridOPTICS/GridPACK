// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   powerflow_solver.hpp
 * @author William A. Perkins
 * @date   2013-04-29 07:23:10 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------

#ifndef _powerflow_solver_hpp_
#define _powerflow_solver_hpp_

#include "gridpack/parallel/distributed.hpp"
#include "gridpack/utility/uncopyable.hpp"
#include "gridpack/utility/matrix.hpp"

namespace gridpack {

// -------------------------------------------------------------
//  class PowerflowSolver
// -------------------------------------------------------------
/// 
/**
 * This abstract class encapsulates the solution of the power flow
 * problem.
 *
 * This class is instantiated using a complete, previously prepared
 * admittance matrix.  
 *
 * Collective on the communicator used by the admittance matrix
 * 
 */
class PowerflowSolver 
  : public parallel::Distributed,
    public utility::Configurable,
    public utility::SanityInterface,
    private utility::Uncopyable
{ 
public:

  /// Default constructor.
  PowerflowSolver(math::Matrix *admittance, 
                  math::Vector *power_injection, 
                  math::Vector *bus_voltage);

  /// Destructor
  ~PowerflowSolver(void);

  /// Solve the power flow equations given power injection and estimated complex voltage
  void solve(void)
  {
    this->solve_();
  }

  /// Fill the specified Vector with the current solution estimate (::bus_voltage_)
  void solution(math::Vector& x) const
  {
    x.copy(*bus_voltage_);
  }

protected:

  /// Admittance matrix (not modified)
  const math::Matrix *admittance_;

  /// Bus power injections (not modified)
  const math::Vector *power_injection_;
  
  /// Complex bus voltage (estimated)
  math::Vector *bus_voltage_;

  /// Performs initial set system setup
  virtual void solve_setup_(void) = 0;

  /// Performs 
  virtual void solve_(void) = 0;

};


} // namespace gridpack


#endif
