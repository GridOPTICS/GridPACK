/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all psssim parameters that can be read in from
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

#ifndef _PSSSIM_VAR_HPP_
#define _PSSSIM_VAR_HPP_

// PSSSIM DATA
/**
 * Flag to indicate that PSSSIM is present
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
 * PSSSIM inputtype
 * type: integer
 * indexed
 */
#define PSSSIM_INPUTTYPE "PSSSIM_INPUTTYPE"

/**
 * PSSSIM bus1 
 * type: ineger
 * indexed
 */
#define PSSSIM_BUS1 "PSSSIM_BUS1"

/**
 * PSSSIM bus2 
 * type: ineger
 * indexed
 */
#define PSSSIM_BUS2 "PSSSIM_BUS2"

/**
 * PSSSIM bus3 
 * type: ineger
 * indexed
 */
#define PSSSIM_BUS3 "PSSSIM_BUS3"

/**
 * PSSSIM bus4 
 * type: ineger
 * indexed
 */
#define PSSSIM_BUS4 "PSSSIM_BUS4"

/**
 * PSSSIM bus5 
 * type: ineger
 * indexed
 */
#define PSSSIM_BUS5 "PSSSIM_BUS5"

/**
 * PSSSIM bus6 
 * type: ineger
 * indexed
 */
#define PSSSIM_BUS6 "PSSSIM_BUS6"

/**
 * PSSSIM gainK 
 * type: real float
 * indexed
 */
#define PSSSIM_GAINK "PSSSIM_GAINK"

/**
 * PSSSIM TW 
 * type: real float
 * indexed
 */
#define PSSSIM_TW "PSSSIM_TW"

/**
 * PSSSIM T1
 * type: real float
 * indexed
 */
#define PSSSIM_T1 "PSSSIM_T1"

/**
 * PSSSIM T2
 * type: real float
 * indexed
 */
#define PSSSIM_T2 "PSSSIM_T2"

/**
 * PSSSIM T3
 * type: real float
 * indexed
 */
#define PSSSIM_T3 "PSSSIM_T3"

/**
 * PSSSIM T4
 * type: real float
 * indexed
 */
#define PSSSIM_T4 "PSSSIM_T4"

/**
 * PSSSIM MAXOUT
 * type: real float
 * indexed
 */
#define PSSSIM_MAXOUT "PSSSIM_MAXOUT"

/**
 * PSSSIM MINOUT
 * type: real float
 * indexed
 */
#define PSSSIM_MINOUT "PSSSIM_MINOUT"

#endif /* _PSSSIM_VAR_HPP_ */
