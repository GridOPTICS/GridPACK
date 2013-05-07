// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   math.hpp
 * @author William A. Perkins
 * @date   2013-05-03 14:13:47 d3g096
 * 
 * @brief  
 * 
 * 
 */
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
