/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all wind machine parameters that can be
 * read in from PTI format files. Each parameter has a corresponding macro
 * that can be used as a unique string to identify the parameter. The use of
 * macros instead of using strings directly will provide extra safety by
 * forcing compiler errors in the case of typos or spelling mistakes.
 */

/**
 *  Variables that can be associated more than once for a bus or a branch can be
 *  indexed by an integer to distinguish different instances. For example,
 *  multiple generators can be associated with a bus and multiple transmission
 *  elements can be associated with a branch. The variables that have an associated
 *  index are denoted with the keyword "indexed".
 */

#ifndef _WIND_VAR_HPP_
#define _WIND_VAR_HPP_

// WIND MODEL DATA
/**
 * Flag to indicate that WIND aerodynamic model is present
 * type: boolean
 * indexed
 */
#define HAS_WIND_AERODYNAMIC "HAS_WIND_AERODYNAMIC"

/**
 * Flag to indicate that WIND drive train model is present
 * type: boolean
 * indexed
 */
#define HAS_WIND_DRIVETRAIN "HAS_WIND_DRIVETRAIN"

/**
 * Flag to indicate that WIND pitch control model is present
 * type: boolean
 * indexed
 */
#define HAS_WIND_PITCHCONTROL "HAS_WIND_PITCHCONTROL"

/**
 * Flag to indicate that WIND torqu control model is present
 * type: boolean
 * indexed
 */
#define HAS_WIND_TORQUECONTROL "HAS_WIND_TORQUECONTROL"

/**
 * Wind aero-dynamic model
 * type: string
 * indexed
 */
#define WIND_AERODYNAMIC "WIND_AERODYNAMIC"

/**
 * Wind drive train model
 * type: string
 * indexed
 */
#define WIND_DRIVETRAIN "WIND_DRIVETRAIN"

/**
 * Wind pitch control model
 * type: string
 * indexed
 */
#define WIND_PITCHCONTROL "WIND_PITCHCONTROL"

/**
 * Wind torque control model
 * type: string
 * indexed
 */
#define WIND_TORQUECONTROL "WIND_TORQUECONTROL"

/**
 * Wind model ID
 * type: string
 * indexed
 */
#define WIND_ID "WIND_ID"

/**
 * Parameters for different components of the wind model are
 * identified by two-character tags
 * Aerodynamic: AD
 * Drive train: DT
 * Pitch control: PC
 * Torque control: TC
 */

/**
 * Wind model ICONM
 * type: integer
 * indexed
 */
#define WIND_TC_TFLAG "WIND_TC_TFLAG"

/**
 * Flag for speed control (0)/power control (1)
 * type: real float
 * indexed
 */
#define WIND_AD_KA "WIND_AD_KA"

/**
 * Initial pitch angle
 * type: real float
 * indexed
 */
#define WIND_AD_THETA "WIND_AD_THETA"

/**
 * Total inertial constant
 * type: real float
 * indexed
 */
#define WIND_DT_H "WIND_DT_H"

/**
 * Machine damping factor
 * type: real float
 * indexed
 */
#define WIND_DT_DAMP "WIND_DT_DAMP"

/**
 * Turbine inertial fraction
 * type: real float
 * indexed
 */
#define WIND_DT_HFRAC "WIND_DT_HFRAC"

/**
 * First shaft torsional resonant frequency
 * type: real float
 * indexed
 */
#define WIND_DT_FREQ1 "WIND_DT_FREQ1"

/**
 * Shaft damping factor
 * type: real float
 * indexed
 */
#define WIND_DT_DSHAFT "WIND_DT_DSHAFT"

/**
 * Pitch-control integral gain
 * type: real float
 * indexed
 */
#define WIND_PC_KIW "WIND_PC_KIW"

/**
 * Pitch-control proportional gain
 * type: real float
 * indexed
 */
#define WIND_PC_KPW "WIND_PC_KPW"

/**
 * Pitch-compensation integral gain
 * type: real float
 * indexed
 */
#define WIND_PC_KIC "WIND_PC_KIC"

/**
 * Pitch-compensation proportional gain
 * type: real float
 * indexed
 */
#define WIND_PC_KPC "WIND_PC_KPC"

/**
 * Gain
 * type: real float
 * indexed
 */
#define WIND_PC_KCC "WIND_PC_KCC"

/**
 * Blade response time constant
 * type: real float
 * indexed
 */
#define WIND_PC_TP "WIND_PC_TP"

/**
 * Maximum pitch angle
 * type: real float
 * indexed
 */
#define WIND_PC_THETAMAX "WIND_PC_THETAMAX"

/**
 * Minimum pitch angle
 * type: real float
 * indexed
 */
#define WIND_PC_THETAMIN "WIND_PC_THETAMIN"

/**
 * Maximum pitch angle rate
 * type: real float
 * indexed
 */
#define WIND_PC_RTHETAMAX "WIND_PC_RTHETAMAX"

/**
 * Minimum pitch angle rate
 * type: real float
 * indexed
 */
#define WIND_PC_RTHETAMIN "WIND_PC_RTHETAMIN"

/**
 * Proportional gain in torque regulator
 * type: real float
 * indexed
 */
#define WIND_TC_KPP "WIND_TC_KPP"

/**
 * Integrator gain in torque regulator
 * type: real float
 * indexed
 */
#define WIND_TC_KIP "WIND_TC_KIP"

/**
 * Electrical power filter time constant
 * type: real float
 * indexed
 */
#define WIND_TC_TP "WIND_TC_TP"

/**
 * Speed reference time constant
 * type: real float
 * indexed
 */
#define WIND_TC_TWREF "WIND_TC_TWREF"

/**
 * Max limit in torque regulator
 * type: real float
 * indexed
 */
#define WIND_TC_TEMAX "WIND_TC_TEMAX"

/**
 * Min limit in torque regulator
 * type: real float
 * indexed
 */
#define WIND_TC_TEMIN "WIND_TC_TEMIN"

/**
 * Power
 * type: real float
 * indexed
 */
#define WIND_TC_P1 "WIND_TC_P1"

/**
 * Shaft speed for power p1
 * type: real float
 * indexed
 */
#define WIND_TC_SPD1 "WIND_TC_SPD1"

/**
 * Power
 * type: real float
 * indexed
 */
#define WIND_TC_P2 "WIND_TC_P2"

/**
 * Shaft speed for power p2
 * type: real float
 * indexed
 */
#define WIND_TC_SPD2 "WIND_TC_SPD2"

/**
 * Power
 * type: real float
 * indexed
 */
#define WIND_TC_P3 "WIND_TC_P3"

/**
 * Shaft speed for power p3
 * type: real float
 * indexed
 */
#define WIND_TC_SPD3 "WIND_TC_SPD3"

/**
 * Power
 * type: real float
 * indexed
 */
#define WIND_TC_P4 "WIND_TC_P4"

/**
 * Shaft speed for power p4
 * type: real float
 * indexed
 */
#define WIND_TC_SPD4 "WIND_TC_SPD4"

/**
 * Total turbine rating
 * type: real float
 * indexed
 */
#define WIND_TC_TRATE "WIND_TC_TRATE"

#endif /* _WIND_VAR_HPP_ */
