/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for shunt parameters that can be read in from
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

#ifndef SHUNT_VAR_HPP_
#define SHUNT_VAR_HPP_

// Shunt data

/**
 * Number of fixed shunts on the bus
 * type: integer
 */
#define SHUNT_NUMBER "SHUNT_NUMBER"

/**
 * Character string to identify shunt on buses with more than one shunt
 * type: string
 * indexed
 */
#define SHUNT_ID "SHUNT_ID"

/**
 * Flag for shunt being on or off
 * type: integer
 * indexed
 */
#define SHUNT_STATUS "SHUNT_STATUS"

/**
 * Active component of the shunt admittance to ground, entered in MW. 
 * Default value is 0.0
 * type: real float
 * indexed
 */
#define BUS_SHUNT_GL "BUS_SHUNT_GL"

/**
 * Reactive component of shunt admittance to ground, entered in Mvar. 
 * positive for a capacitor and negative for a reactor
 * Default value is 0.0
 * type: real float
 * indexed
 */
#define BUS_SHUNT_BL "BUS_SHUNT_BL"

// SWITCHED SHUNT DATA

/**
 * Bus number to which the shunt is connected
 * type: integer
 */
#define SWSHUNT_BUSNUMBER "SWSHUNT_BUSNUMBER"

/**
 * Bus number to which fixed shunt is connected
 * type: integer
 */
#define SHUNT_BUSNUMBER "SHUNT_BUSNUMBER"


/**
 * Control mode:
 *  0 - fixed
 *  1 - discrete adjustment, controlling voltage locally or at bus SHUNT_SWREM
 *  2 - continuous adjustment, controlling voltage locally or at bus SHUNT_SWREM
 *  3 - discrete adjustment, controlling reactive power output of the plant at bus SHUNT_SWREM
 *  4 - discrete adjustment, controlling reactive power output of the VSC dc line converter at bus SHUNT_SWREM of the VSC dc line whose name is specified as SHUNT_RMIDNT
 *  5 - discrete adjustment, controlling admittance setting of the switched shunt at bus SHUNT_SWREM.
 * SHUNT_MODSW = 1 by default.
 * type: integer
 */
#define SHUNT_MODSW "SHUNT_MODSW"

/**
 * Adjustment method:
 *   0 - steps and blocks are switched on in iput order, and off in reverse
 *   order
 *   1 - steps and blocks are switched on and off such that the next highest (or
 *   lowest, as appropriate) total admittance is achieved.
 *   SHUNT_ADJM = 0 by default
 *   type: integer
 */
#define SHUNT_ADJM "SHUNT_ADJM"

/**
 * Initial switched shunt status of one for in-service and zero for
 * out-of-service
 * SHUNT_SWCH_STAT = 1 by default
 * type: integer
 */
#define SHUNT_SWCH_STAT "SHUNT_SWCH_STAT"

/**
 * When SHUNT_MODSW is 1 or 2, the controlled voltage upper limit; entered in pu.
 * When SHUNT_MODSW is 3, 4 or 5, the controlled reactive power range upper limit;
 * entered in pu of the total reactive power range of the controlled voltage
 * controlling device.
 * SHUNT_VSWHI = 1.0 by default.
 * type: real float
 */
#define SHUNT_VSWHI "SHUNT_VSWHI"

/**
 * When SHUNT_MODSW is 1 or 2, the controlled voltage lower limit; entered in pu.
 * When SHUNT_MODSW is 3, 4 or 5, the controlled reactive power range lower limit;
 * entered in pu of the total reactive power range of the controlled voltage
 * controlling device. SHUNT_VSWLO = 1.0 by default
 * type: real float
 */
#define SHUNT_VSWLO "SHUNT_VSWLO"

/**
 * Bus number of the bus whose voltage or connected equipment reactive power output
 * is controlled by this switched shunt
 * type: integer
 */
#define SHUNT_SWREM "SHUNT_SWREM"
#define SHUNT_SWREG "SHUNT_SWREG"

/**
 * A node number of bus SWREG
 * type: integer
 */
#define SHUNT_NREG "SHUNT_NREG"

/**
 * Percent of the total MVar required to hold the voltage at the bus controlled
 * by bus SWSHUNT_BUSNUMBER that are contributed by this switched shunt;
 * SHUNT_RMPCT must be positive.
 * SHUNT_RMPCT = 100.0 by default
 * type: float
 */
#define SHUNT_RMPCT "SHUNT_RMPCT"

/**
 * When SHUNT_MODSW is 4, the name of the VSC dc line where the converter bus is
 * specified in SHUNT_SWREM.
 * When SHUNT_MODSW is 6, the name of the FACTS device where the shunt element's
 * reactive output is to be controlled.
 * SHUNT_RMIDNT is blank by default
 * type: string
 */
#define SHUNT_RMIDNT "SHUNT_RMIDNT"

/**
 * Initial switched shunt admittance; entered in Mvar at unity voltage
 * Default value: 0.0
 * type: real float
 */
#define SHUNT_BINIT "SHUNT_BINIT"

/**
 * Initial switched shunt status for block i
 * Default value: 1
 * type: integer
 */
#define SHUNT_S1 "SHUNT_S1"
#define SHUNT_S2 "SHUNT_S2"
#define SHUNT_S3 "SHUNT_S3"
#define SHUNT_S4 "SHUNT_S4"
#define SHUNT_S5 "SHUNT_S5"
#define SHUNT_S6 "SHUNT_S6"
#define SHUNT_S7 "SHUNT_S7"
#define SHUNT_S8 "SHUNT_S8"

/**
 * Number of steps for block i. The first zero value of Ni or Bi is interpreted as
 * the end of the switched shunt blocks for bus I.
 * Default value: 0
 * type: integer
 */
#define SHUNT_N1 "SHUNT_N1"
#define SHUNT_N2 "SHUNT_N2"
#define SHUNT_N3 "SHUNT_N3"
#define SHUNT_N4 "SHUNT_N4"
#define SHUNT_N5 "SHUNT_N5"
#define SHUNT_N6 "SHUNT_N6"
#define SHUNT_N7 "SHUNT_N7"
#define SHUNT_N8 "SHUNT_N8"

/**
 * Admittance increment for each of Ni steps in block i; entered in Mvar at unity voltage.
 * Default value: 0
 * type: real float
 */
#define SHUNT_B1 "SHUNT_B1"
#define SHUNT_B2 "SHUNT_B2"
#define SHUNT_B3 "SHUNT_B3"
#define SHUNT_B4 "SHUNT_B4"
#define SHUNT_B5 "SHUNT_B5"
#define SHUNT_B6 "SHUNT_B6"
#define SHUNT_B7 "SHUNT_B7"
#define SHUNT_B8 "SHUNT_B8"

#endif /* _SHUNT_VAR_HPP_ */
