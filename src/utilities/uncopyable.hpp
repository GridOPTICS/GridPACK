//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   uncopyable.hpp
 * @author William A. Perkins
 * @date   2013-05-07 15:12:17 d3g096
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
// Created April  1, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _uncopyable_hpp_
#define _uncopyable_hpp_

#include <boost/utility.hpp>

namespace gridpack {
namespace utility {

typedef boost::noncopyable Uncopyable;

} // namespace math
} // namespace gridpack





#ifdef _DO_INLINE_
#include "uncopyable-inline.hpp"
#endif

#endif
