// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   linear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   Mon Apr  1 08:39:20 2013
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  1, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _linear_solver_implementation_hpp_
#define _linear_solver_implementation_hpp_

#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/utility/noncopyable.hpp"
#include "gridpack/utility/configurable.hpp"
#include "gridpack/parallel/distributed.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearSolverImplementation
// -------------------------------------------------------------
class LinearSolverImplementation 
  : public parallel::Distributed,
    public parallel::UnCopyable,
    public utility::Configurable
{
public:
  
  /// Default constructor.
  LinearSolverImplementation(const parallel::Distribution& dist,
                             const Matrix& A);
  
  /// Destructor
  ~LinearSolverImplementation(void);
  
  /// Solve w/ the specified RHS, put result in specified vector
  void solve(const Vector& b, Vector& x) const
  {
    this->solve_(b, x);
  }

  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor)
  {
    this->accept_(visitor);
  }


protected:

  /// The coefficient matrix (may not need to remember)
  const Matrix& a_;
  
  /// Solve w/ the specified RHS, put result in specified vector
  /** 
   * Can be called repeatedly with different @c b and @c x vectors
   * 
   * @param b Vector containing right hand side of linear system
   * @param x Vector containing initial solution estimate; final
   * solution is put into this.
   */
  virtual void solve_(const Vector& b, Vector& x) const = 0;

  /// Allow visits by implementation visitors
  virtual void accept_(ImplementationVisitor& visitor) = 0;

};

} // namespace math
} // namespace gridpack



#endif
