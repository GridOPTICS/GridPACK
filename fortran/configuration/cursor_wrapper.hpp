// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   cursor_wrapper.hpp
 * @author William A. Perkins
 * @date   2014-08-14 14:51:10 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 14, 2014 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _cursor_wrapper_hpp_
#define _cursor_wrapper_hpp_

#include <gridpack/configuration/configuration.hpp>

// -------------------------------------------------------------
// CursorWrapper
// -------------------------------------------------------------
struct CursorWrapper {
  gridpack::utility::Configuration::CursorPtr cursor;
};

#endif
