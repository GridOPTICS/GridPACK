// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dae_solver_functions.hpp
 * @author William A. Perkins
 * @date   2013-11-11 09:54:39 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#ifndef _dae_solver_functions_hpp_
#define _dae_solver_functions_hpp_

#include <boost/function.hpp>
#include <gridpack/math/matrix.hpp>

namespace gridpack {
namespace math {

typedef 
boost::function<void (const double& time, 
                      const Vector& x, const Vector& xdot, 
                      const double& shift, Matrix& J)> 
DAEJacobianBuilder;

typedef 
boost::function<void (const double& time, 
                      const Vector& x, const Vector& xdot, 
                      Vector& F)> 
DAEFunctionBuilder;



} // namespace math
} // namespace gridpack

#endif
