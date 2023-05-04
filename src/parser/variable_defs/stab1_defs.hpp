/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all stab1 parameters that can be read in from
 * PTI format files. Each parameter has a corresponding macro that can be used
 * as a unique string to identify the parameter. The use of macros instead of
 * using strings directly will provide extra safety by forcing compiler errors
 * in the case of typos or spelling mistakes.
 */

/**
 *  Variables that can be associated more than once for a bus or a branch can be
 *  indexed by an integer to distinguish different instances. For example,
 *  multiple generators can be associated with a bus and multiple transmission
 *  elements can be associated with a branch. The variables that have an associated
 *  index are denoted with the keyword "indexed".
 */

#ifndef _STAB1_VAR_HPP_
#define _STAB1_VAR_HPP_

// STAB1 DATA
/**
 * Flag to indicate that STAB1 is present
 * type: boolean
 * indexed
 */
#define HAS_PSS "HAS_PSS"

/**
 * Exciter model
 * type: string
 * indexed
 */
#define PSS_MODEL "PSS_MODEL"

/**
 * STAB1 J
 * type: real float
 * indexed
 */
#define STAB1_J "STAB1_J"

/**
 * STAB1 J+1
 * type: real float
 * indexed
 */
#define STAB1_J1 "STAB1_J1"

/**
 * STAB1 J+2
 * type: real float
 * indexed
 */
#define STAB1_J2 "STAB1_J2"

/**
 * STAB1 J+3
 * type: real float
 * indexed
 */
#define STAB1_J3 "STAB1_J3"

/**
 * STAB1 J+4
 * type: real float
 * indexed
 */
#define STAB1_J4 "STAB1_J4"

/**
 * STAB1 J+5
 * type: real float
 * indexed
 */
#define STAB1_J5 "STAB1_J5"

/**
 * STAB1 J+6
 * type: real float
 * indexed
 */
#define STAB1_J6 "STAB1_J6"

#endif /* _STAB1_VAR_HPP_ */
