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
 * @date   2015-05-05 10:08:39 d3g096
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

template <typename T, typename I = int>
struct DAEBuilder {
  
  typedef VectorT<T, I> VectorType;
  typedef MatrixT<T, I> MatrixType;

  typedef 
  boost::function<void (const double& time, 
                        const VectorType& x, const VectorType& xdot, 
                        const double& shift, MatrixType& J)> 
  Jacobian;

  typedef 
  boost::function<void (const double& time, 
                        const VectorType& x, const VectorType& xdot, 
                        VectorType& F)> 
  Function;

};


} // namespace math
} // namespace gridpack

#endif
