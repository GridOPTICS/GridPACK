// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_functions.hpp
 * @author William A. Perkins
 * @date   2013-08-12 09:31:54 d3g096
 * 
 * @brief Declaration of function (objects) used by nonlinear solvers
 * to build a Jacobian and RHS.
 * 
 * 
 */
#ifndef _nonlinear_solver_functions_hpp_
#define _nonlinear_solver_functions_hpp_

#include <boost/function.hpp>
#include <gridpack/math/matrix.hpp>

namespace gridpack {
namespace math {

/// A function that is used to build a Jacobian matrix
typedef boost::function<void (const Vector& x, Matrix& J)> JacobianBuilder;

/// A function than is used to build a RHS vector
typedef boost::function<void (const Vector& x, Vector& F)> FunctionBuilder;

} // namespace math
} // namespace gridpack

#endif
