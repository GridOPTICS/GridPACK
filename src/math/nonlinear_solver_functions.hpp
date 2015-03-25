// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_functions.hpp
 * @author William A. Perkins
 * @date   2015-03-25 12:07:07 d3g096
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

// -------------------------------------------------------------
// struct NLSBuilder
// -------------------------------------------------------------
template <typename T, typename I = int>
struct NLSBuilder {

  /// A function object that builds the Jacobian for NonlinearSolver
  typedef boost::function<void (const VectorT<T, I>& x, MatrixT<T, I>& J)> Jacobian;

  /// A function object that builds the RHS for NonlinearSolver
  typedef boost::function<void (const VectorT<T, I>& x, VectorT<T, I>& F)> Function;

};

typedef typename NLSBuilder<ComplexType>::Jacobian ComplexJacobianBuilder;
typedef typename NLSBuilder<ComplexType>::Function ComplexFunctionBuilder;
typedef typename NLSBuilder<RealType>::Jacobian RealJacobianBuilder;
typedef typename NLSBuilder<RealType>::Function RealFunctionBuilder;

/// A function used to build a Jacobian matrix for a nonlinear system of equations
/**
 * This type is used to supply a way to build a Jacobian Matrix to
 * NonlinearSolver instances. The NonlinearSolver will use this each
 * iteration to construct the Matrix @c J from the current solution
 * estimate in @c x. A Jacobian builder may be a function or function
 * object.  For example, this
 *
 @code{.cpp}
 void 
 my_jacobian_builder(const Vector& x, Matrix& J)
 {
 // ...
 J.ready();
 }
 JacobianBuilder j = my_jacobian_builder;
 @endcode
 *
 * and 
 * 
 \code{.cpp}
 struct my_jacobian_builder_type {
 void operator()(const Vector& x, Matrix& J)
 {
 // ...
 J.ready();
 }
 } 
 my_jacobian_builder_type my_jacobian_builder;
 JacobianBuilder j = my_jacobian_builder;
 \endcode
 * 
 * are equivalent. The latter is usually more convenient because
 * necessary information can be included in the function object.
 */
typedef ComplexJacobianBuilder JacobianBuilder;

/// A function used to build the RHS of a nonlinear system of equations
/**
 * This type is used to supply a way to build a right hand side Vector
 * to NonlinearSolver instances. The NonlinearSolver will use this
 * each iteration to construct the Vector @c F from the current
 * solution estimate in @c x. A FunctionBuilder may be a function or
 * function object.  For example, this
 *
 @code{.cpp}
 void 
 my_function_builder(const Vector& x, Vector& F)
 {
 // ...      
 F.ready();
 }
 JacobianBuilder j = my_function_builder;
 @endcode
 *
 * and 
 * 
 \code{.cpp}
 struct my_function_builder_type {
 void operator()(const Vector& x, Vector& F)
 {
 // ...
 F.ready();
 }
 }; 
 my_function_builder_type my_function_builder;
 JacobianBuilder j = boost::ref(my_function_builder);
 \endcode
 * 
 * are equivalent. The latter is usually more convenient because
 * necessary information can be included in the function object.
 */
typedef ComplexFunctionBuilder FunctionBuilder;

} // namespace math
} // namespace gridpack

#endif
