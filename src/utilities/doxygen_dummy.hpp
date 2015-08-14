//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   doxygen_dummy.hpp
 * @author William A. Perkins
 * @date   2015-08-07 09:03:18 d3g096
 * 
 * @brief  Some bogus code tricks to make Doxygen documentation better
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created September 10, 2013 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _doxygen_dummy_hpp_
#define _doxygen_dummy_hpp_

namespace boost { 

/// Bogus class class for boost::shared_ptr
template<class T> class shared_ptr { T *dummy; }; 

/// Bogus class class for boost::scoped_ptr
template<class T> class scoped_ptr { T *dummy; }; 
}

#endif
