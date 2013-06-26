/*
 * This file contains definitions for all parameters that can be read in from
 * PTI format files. Each parameter has a corresponding macro that can be used
 * as a unique string to identify the parameter. The use of macros instead of
 * using strings directly will provide extra safety by forcing compiler errors
 * in the case of typos or spelling mistackes.
 */

#ifndef DICTIONARY_HPP_
#define DICTIONARY_HPP_

/**
 * 0: base case
 * 1: add information to existing case
 * type: integer
 */
#define CASE_ID "CASE_ID"

/**
 * System base MVS. Default value is 100.0
 * type: real float
 */
#define CASE_SBASE "CASE_SBASE"

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
 * Bus base voltage, entered in kV. Default value is 0.0
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
 * Active component of shunt admittance to ground, entered in MW at one per unit
 * voltage. BUS_SHUNT_GL should not include any resistance load, which is
 * entered as part of load data. Default value is 0.0
 * type: real float
 */
#define BUS_SHUNT_GL "BUS_SHUNT_GL"

/**
 * Reactive component of shunt admittance to ground, entered in Mvar at one per
 * unit voltage. BUS_SHUNT_BL should not include any reactive impedence load,
 * which is entered as part of load data, Line charging and line connect shunts
 * which are entered as part of non-transformer branch data, or transformer
 * admittance, which is entered as part of transformer data. BUS_SHUNT_BL is
 * positive for a capacitor and negative for a reactor or an inductive load.
 * Default value is 0.0
 * type: real float
 */
#define BUS_SHUNT_GL "BUS_SHUNT_GL"

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
 * Bus voltage angle, in degrees
 * type: real float
 */
#define BUS_VOLTAGE_ANG "BUS_VOLTAGE_ANG"

/**
 * Owner number
 * type: integer
 */
#define BUS_OWNER "BUS_OWNER"

#endif /* DICTIONARY_HPP_ */
