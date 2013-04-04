// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   math_type.hpp
 * @author William A. Perkins
 * @date   Wed Apr  3 12:12:49 2013
 * 
 * @brief 
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  3, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _math_type_hpp_
#define _math_type_hpp_

#include <complex>

namespace gridpack {
namespace math {

/// The type that all Matrix and Vector instances contain
/**
 * This probably depends on the math library used.  This is OK for
 * PETSc.
 * 
 */
typedef std::complex<double> complex_type;

// If libraries other than PETSc are used, a different complex type
// may be required, along with associated handling methods.

} // namespace utility
} // namespace gridpack

#endif
