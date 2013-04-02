// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   uncopyable.hpp
 * @author William A. Perkins
 * @date   Mon Apr  1 10:09:49 2013
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

typedef boost::noncopyable UnCopyable;

} // namespace math
} // namespace gridpack





#ifdef _DO_INLINE_
#include "uncopyable-inline.hpp"
#endif

#endif
