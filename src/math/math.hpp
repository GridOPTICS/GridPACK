// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   math.hpp
 * @author William A. Perkins
 * @date   2013-04-26 12:07:31 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 26, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _math_hpp_
#define _math_hpp_

namespace gridpack {
namespace math {

#include "gripack/math/math_type.hpp"

/// Do whatever is necessary to initialize the math library
extern void Initialize(void);

/// Do whatever is necessary to shut down the math library
extern void Finalize(void);

} // namespace math
} // namespace gridpack


#endif
