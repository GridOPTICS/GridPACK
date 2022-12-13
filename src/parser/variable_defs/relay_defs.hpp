/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all relay parameters that can be read in from
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

#ifndef _RELAY_VAR_HPP_
#define _RELAY_VAR_HPP_

// RELAY DATA
/**
 * Relay number
 * type: integer
 */
#define RELAY_NUMBER "RELAY_NUMBER"

/**
 * Relay model
 * type: string
 * indexed
 */
#define RELAY_MODEL "RELAY_MODEL"

/**
 * Relay lid
 * type: string
 * indexed
 */
#define RELAY_LID "RELAY_LID"

/**
 * Relay jbus
 * type: integer
 * indexed
 */
#define RELAY_JBUS "RELAY_JBUS"

/**
 * Relay v1
 * type: real float
 * indexed
 */
#define RELAY_V1 "RELAY_V1"

/**
 * Relay t1
 * type: real float
 * indexed
 */
#define RELAY_T1 "RELAY_T1"

/**
 * Relay f1
 * type: real float
 * indexed
 */
#define RELAY_F1 "RELAY_F1"

/**
 * Relay v2
 * type: real float
 * indexed
 */
#define RELAY_V2 "RELAY_V2"

/**
 * Relay t2
 * type: real float
 * indexed
 */
#define RELAY_T2 "RELAY_T2"

/**
 * Relay f2
 * type: real float
 * indexed
 */
#define RELAY_F2 "RELAY_F2"

/**
 * Relay v3
 * type: real float
 * indexed
 */
#define RELAY_V3 "RELAY_V3"

/**
 * Relay t3
 * type: real float
 * indexed
 */
#define RELAY_T3 "RELAY_T3"

/**
 * Relay f3
 * type: real float
 * indexed
 */
#define RELAY_F3 "RELAY_F3"

/**
 * Relay tb
 * type: real float
 * indexed
 */
#define RELAY_TB "RELAY_TB"

/**
 * Relay mins
 * type: integer
 * indexed
 */
#define RELAY_MINS "RELAY_MINS"

/**
 * Relay frequency bus
 * type: integer
 * indexed
 */
#define RELAY_FREBUS "RELAY_FREBUS"

/**
 * Generator ID
 * type: string
 * indexed
 */
#define RELAY_GENID "RELAY_GENID"

/**
 * Relay fl
 * type: real float
 * indexed
 */
#define RELAY_FL "RELAY_FL"

/**
 * Relay fu
 * type: real float
 * indexed
 */
#define RELAY_FU "RELAY_FU"

/**
 * Relay tp
 * type: real float
 * indexed
 */
#define RELAY_TP "RELAY_TP"

/**
 * Relay ID
 * type: string
 * indexed
 */
#define RELAY_ID "RELAY_ID"

/**
 * Relay rs
 * type: integer
 * indexed
 */
#define RELAY_RS "RELAY_RS"

/**
 * Relay mtype
 * type: integer
 * indexed
 */
#define RELAY_MTYPE "RELAY_MTYPE"

/**
 * Relay bmon
 * type: integer
 * indexed
 */
#define RELAY_BMON "RELAY_BMON"

/**
 * Relay ibus1
 * type: integer
 * indexed
 */
#define RELAY_IBUS1 "RELAY_IBUS1"

/**
 * Relay jbus1
 * type: integer
 * indexed
 */
#define RELAY_JBUS1 "RELAY_JBUS1"

/**
 * Relay ID1
 * type: string
 * indexed
 */
#define RELAY_ID1 "RELAY_ID1"

/**
 * Relay ibus2
 * type: integer
 * indexed
 */
#define RELAY_IBUS2 "RELAY_IBUS2"

/**
 * Relay jbus2
 * type: integer
 * indexed
 */
#define RELAY_JBUS2 "RELAY_JBUS2"

/**
 * Relay ID2
 * type: string
 * indexed
 */
#define RELAY_ID2 "RELAY_ID2"

/**
 * Relay ibus3
 * type: integer
 * indexed
 */
#define RELAY_IBUS3 "RELAY_IBUS3"

/**
 * Relay jbus3
 * type: integer
 * indexed
 */
#define RELAY_JBUS3 "RELAY_JBUS3"

/**
 * Relay ID2
 * type: string
 * indexed
 */
#define RELAY_ID3 "RELAY_ID3"

/**
 * Relay Zone1 time
 * type: real float
 * indexed
 */
#define RELAY_ZONE1_TIME "RELAY_ZONE1_TIME"

/**
 * Relay Zone1 reach
 * type: real float
 * indexed
 */
#define RELAY_ZONE1_REACH "RELAY_ZONE1_REACH"

/**
 * Relay Zone1 cenang
 * type: real float
 * indexed
 */
#define RELAY_ZONE1_CENANG "RELAY_ZONE1_CENANG"

/**
 * Relay Zone1 cenang
 * type: real float
 * indexed
 */
#define RELAY_ZONE1_CENDIS "RELAY_ZONE1_CENDIS"

/**
 * Relay Zone2 time
 * type: real float
 * indexed
 */
#define RELAY_ZONE2_TIME "RELAY_ZONE2_TIME"

/**
 * Relay Zone2 reach
 * type: real float
 * indexed
 */
#define RELAY_ZONE2_REACH "RELAY_ZONE2_REACH"

/**
 * Relay Zone2 cenang
 * type: real float
 * indexed
 */
#define RELAY_ZONE2_CENANG "RELAY_ZONE2_CENANG"

/**
 * Relay Zone2 cenang
 * type: real float
 * indexed
 */
#define RELAY_ZONE2_CENDIS "RELAY_ZONE2_CENDIS"

/**
 * Relay Zone3 time
 * type: real float
 * indexed
 */
#define RELAY_ZONE3_TIME "RELAY_ZONE3_TIME"

/**
 * Relay Zone3 reach
 * type: real float
 * indexed
 */
#define RELAY_ZONE3_REACH "RELAY_ZONE3_REACH"

/**
 * Relay Zone2 cenang
 * type: real float
 * indexed
 */
#define RELAY_ZONE3_CENANG "RELAY_ZONE3_CENANG"

/**
 * Relay Zone3 cenang
 * type: real float
 * indexed
 */
#define RELAY_ZONE3_CENDIS "RELAY_ZONE3_CENDIS"

/**
 * Relay dirang
 * type: real float
 * indexed
 */
#define RELAY_DIRANG "RELAY_DIRANG"

/**
 * Relay thcur
 * type: real float
 * indexed
 */
#define RELAY_THCUR "RELAY_THCUR"

/**
 * Relay sebtime
 * type: real float
 * indexed
 */
#define RELAY_SEBTIME "RELAY_SEBTIME"

/**
 * Relay serctime
 * type: real float
 * indexed
 */
#define RELAY_SERCTIME "RELAY_SERCTIME"

/**
 * Relay trbtime
 * type: real float
 * indexed
 */
#define RELAY_TRBTIME "RELAY_TRBTIME"

/**
 * Relay trrctime
 * type: real float
 * indexed
 */
#define RELAY_TRRCTIME "RELAY_TRRCTIME"

/**
 * Relay bltype1
 * type: integer
 * indexed
 */
#define RELAY_BLTYPE1 "RELAY_BLTYPE1"

/**
 * Relay blint1
 * type: real float
 * indexed
 */
#define RELAY_BLINT1 "RELAY_BLINT1"

/**
 * Relay blro1
 * type: real float
 * indexed
 */
#define RELAY_BLRO1 "RELAY_BLRO1"

/**
 * Relay bltype2
 * type: integer
 * indexed
 */
#define RELAY_BLTYPE2 "RELAY_BLTYPE2"

/**
 * Relay blint2
 * type: real float
 * indexed
 */
#define RELAY_BLINT2 "RELAY_BLINT2"

/**
 * Relay blro2
 * type: real float
 * indexed
 */
#define RELAY_BLRO2 "RELAY_BLRO2"

#endif /* DICTIONARY_HPP_ */
