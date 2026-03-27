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

// IEEEST PSS DATA
#define IEEEST_MODE  "IEEEST_MODE"
#define IEEEST_A1    "IEEEST_A1"
#define IEEEST_A2    "IEEEST_A2"
#define IEEEST_A3    "IEEEST_A3"
#define IEEEST_A4    "IEEEST_A4"
#define IEEEST_A5    "IEEEST_A5"
#define IEEEST_A6    "IEEEST_A6"
#define IEEEST_T1    "IEEEST_T1"
#define IEEEST_T2    "IEEEST_T2"
#define IEEEST_T3    "IEEEST_T3"
#define IEEEST_T4    "IEEEST_T4"
#define IEEEST_T5    "IEEEST_T5"
#define IEEEST_T6    "IEEEST_T6"
#define IEEEST_KS    "IEEEST_KS"
#define IEEEST_LSMAX "IEEEST_LSMAX"
#define IEEEST_LSMIN "IEEEST_LSMIN"
#define IEEEST_VCU   "IEEEST_VCU"
#define IEEEST_VCL   "IEEEST_VCL"

// ST2CUT PSS DATA
#define ST2CUT_MODE   "ST2CUT_MODE"
#define ST2CUT_MODE2  "ST2CUT_MODE2"
#define ST2CUT_K1     "ST2CUT_K1"
#define ST2CUT_K2     "ST2CUT_K2"
#define ST2CUT_T1     "ST2CUT_T1"
#define ST2CUT_T2     "ST2CUT_T2"
#define ST2CUT_T3     "ST2CUT_T3"
#define ST2CUT_T4     "ST2CUT_T4"
#define ST2CUT_T5     "ST2CUT_T5"
#define ST2CUT_T6     "ST2CUT_T6"
#define ST2CUT_T7     "ST2CUT_T7"
#define ST2CUT_T8     "ST2CUT_T8"
#define ST2CUT_T9     "ST2CUT_T9"
#define ST2CUT_T10    "ST2CUT_T10"
#define ST2CUT_LSMAX  "ST2CUT_LSMAX"
#define ST2CUT_LSMIN  "ST2CUT_LSMIN"
#define ST2CUT_VCU    "ST2CUT_VCU"
#define ST2CUT_VCL    "ST2CUT_VCL"

#endif /* _PSSSIM_VAR_HPP_ */
