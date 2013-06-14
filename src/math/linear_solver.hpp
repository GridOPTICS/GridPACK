// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   linear_solver.hpp
 * @author William A. Perkins
 * @date   2013-06-14 12:04:59 d3g096
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

// SCCS ID: $Id$ Battelle PNL

#ifndef _linear_solver_hpp_
#define _linear_solver_hpp_

#include <boost/scoped_ptr.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/math/linear_solver_implementation.hpp>

namespace gridpack {
namespace math {

class ImplementationVisitor;

// -------------------------------------------------------------
//  class LinearSolver
// -------------------------------------------------------------
class LinearSolver 
  : public parallel::Distributed,
    private utility::Uncopyable
{
public:
  
  /// Default constructor.
  LinearSolver(const Matrix& A);
  
  /// Destructor
  ~LinearSolver(void);

  /// Solve w/ the specified RHS, put result in specified vector
  void solve(const Vector& b, Vector& x) const
  {
    p_solver->solve(b, x);
  }

  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor)
  {
    p_solver->accept(visitor);
  }

  /// Allow visits by implemetation visitor
  void accept(ConstImplementationVisitor& visitor) const
  {
    p_solver->accept(visitor);
  }

protected:

  /// Where the work really happens
  boost::scoped_ptr<LinearSolverImplementation> p_solver;

};

} // namespace math
} // namespace gridpack

#endif
