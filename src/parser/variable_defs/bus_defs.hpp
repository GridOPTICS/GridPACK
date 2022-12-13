/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all bus parameters that can be read in from
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

#ifndef _BUS_VAR_HPP_
#define _BUS_VAR_HPP_

// BUS DATA

/**
 * Bus number (1 though 999999)
 * type: integer
 */
#define BUS_NUMBER "BUS_NUMBER"

/**
 * Alpha-numeric identifier assigned to bus
 * type: string
 */
#define BUS_NAME "BUS_NAME"

/**
 * Alpha-numeric identifier assigned to new buses that may be
 * as part of the parsing process
 * type: string
 */
#define NEW_BUS_TYPE "NEW_BUS_TYPE"

/**
 * Bus base voltage, entered in kV. 
 * type: real float
 */
#define BUS_BASEKV "BUS_BASEKV"

/**
 * Bus type
 * 1: load bus
 * 2: generator bus
 * 3: swing bus
 * 4: isolated bus
 * type: integer
 */
#define BUS_TYPE "BUS_TYPE"

/**
 * Area number
 * type: integer
 */
#define BUS_AREA "BUS_AREA"

/**
 * Zone number
 * type: integer
 */
#define BUS_ZONE "BUS_ZONE"

/**
 * Bus voltage magnitude, in p.u.
 * type: real float
 */
#define BUS_VOLTAGE_MAG "BUS_VOLTAGE_MAG"

/**
 * Bus voltage phase angle, in degrees
 * type: real float
 */
#define BUS_VOLTAGE_ANG "BUS_VOLTAGE_ANG"

/**
 * Maximum allowable bus voltage magnitude, in p. u.
 * type: real float
 */
#define BUS_VOLTAGE_MAX "BUS_VOLTAGE_MAX"

/**
 * Minimum allowable bus voltage magnitude, in p. u.
 * type: real float
 */
#define BUS_VOLTAGE_MIN "BUS_VOLTAGE_MIN"

/**
 * Owner number
 * type: integer
 */
#define BUS_OWNER "BUS_OWNER"

/**
 * Flag that indicates that bus was generated from a 3-winding transformer
 * type: boolean
 */
#define BUS_3WINDING "BUS_3WINDING"

// Bus sequence parameters
/**
 * Active component of negative sequance shunt addmittance
 * type: real float
 */
#define BUS_SEQ_GNEG "BUS_SEQ_GNEG"

/**
 * Reactive component of negative sequance shunt addmittance
 * type: real float
 */
#define BUS_SEQ_BNEG "BUS_SEQ_BNEG"

/**
 * Active component of zero sequance shunt addmittance
 * type: real float
 */
#define BUS_SEQ_GZERO "BUS_SEQ_GZERO"

/**
 * Reactive component of zero sequance shunt addmittance
 * type: real float
 */
#define BUS_SEQ_BZERO "BUS_SEQ_BZERO"

#endif /* _BUS_VAR_HPP_ */
