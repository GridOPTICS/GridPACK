// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   uncopyable.hpp
 * @author William A. Perkins
 * @date   2013-10-09 13:44:11 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

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
