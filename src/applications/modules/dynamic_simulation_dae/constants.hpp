/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsimnetwork.cpp
 * @author Shrirang Abhyankar
 * @date   02/17/19
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _constants_h
#define _constants_h

#include <gridpack/include/gridpack.hpp>

#define PI (4.0*atan(1.0))
#define FREQ 60.0
#define OMEGA_S (2.0*PI*60.0)
#define DEFAULT_MVABASE 100.0 // Default MVA base 100 MVA

enum DSMode{NONE,INIT_X,RESIDUAL_EVAL,XVECTOBUS,XDOTVECTOBUS,FAULT_EVAL,XVECPRETOBUS};

#endif
