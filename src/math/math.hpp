// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   math.hpp
 * @author William A. Perkins
 * @date   2013-05-08 15:10:43 d3g096
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

// Do whatever is necessary to initialize the math library
extern void Initialize(void);

/// Is the math library initialized?
extern bool Initialized(void);

/// Do whatever is necessary to shut down the math library
extern void Finalize(void);

} // namespace math
} // namespace gridpack


#endif
