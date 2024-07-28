/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   constants.hpp
 * 
 * @brief  Declares all the constants
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _constants_h
#define _constants_h

#include <gridpack/include/gridpack.hpp>

#define PI 3.1415926535897932384626433832795029
#define FREQ 60.0
#define OMEGA_S (2.0*PI*60.0)
#define DEFAULT_MVABASE 100.0 // Default MVA base 100 MVA
#define TWO_OVER_THREE (2.0/3.0)
#define TWOPI_OVER_THREE (TWO_OVER_THREE*PI)

// This is used during calculations to set x/xdot vector to bus, calculate residual, 
enum EMTMode{NONE,INIT_X,RESIDUAL_EVAL,XVECTOBUS,XDOTVECTOBUS,FAULT_EVAL,XVECPRETOBUS};

// Type of integration algorithm used for machines
// Explicit - Use explicit integration scheme (forward Euler)
// Implicit - Use implicit integration scheme (machine equations included in the system equations)
// IMPLICITEXPLICIT - Uses PETSc's IMEX method, machines use explicit, network uses implicit
// Note :- Some of the machines may not have all the different variants available
enum EMTMachineIntegrationType{EXPLICIT,IMPLICIT,IMPLICITEXPLICIT};

// Buffer size for exchange.
// Note: GridPACK does not allow different buffer sizes for each
// bus and branch components. So we need to hack it to use equal sizes. Hence, using this rather large buffer size for buses and branches
#define BRANCHBUFSIZE 60

#endif
