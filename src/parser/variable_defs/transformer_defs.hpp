/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all transformer parameters that can be
 * read in from PTI format files. Each parameter has a corresponding macro that
 * can be used as a unique string to identify the parameter. The use of macros
 * instead of using strings directly will provide extra safety by forcing
 * compiler errors in the case of typos or spelling mistakes.
 */

/**
 *  Variables that can be associated more than once for a bus or a branch can be
 *  indexed by an integer to distinguish different instances. For example,
 *  multiple generators can be associated with a bus and multiple transmission
 *  elements can be associated with a branch. The variables that have an associated
 *  index are denoted with the keyword "indexed".
 */

#ifndef _TRANSFORMER_VAR_HPP_
#define _TRANSFORMER_VAR_HPP_

// TRANSFORMER DATA 
/**
 * Not parsed in PTI v23
 * Bus number to which the first winding is connected.
 * type: integer
 */
#define TRANSFORMER_BUS1 "TRANSFORMER_BUS1"

/**
 * Not parsed in PTI v23
 * Bus number to which the second winding is connected
 * type: integer
 */
#define TRANSFORMER_BUS2 "TRANSFORMER_BUS2"

/**
 * Not parsed in PTI v23
 * Non-blank alphanumeric transformer circuit identifier
 * type: string
 */
#define TRANSFORMER_CKT "TRANSFORMER_CKT"

/**
 * Number of bus to control. If different from BUS1 or BUS2 then sign determines
 * control. Positive sign, close to impedance (untapped) bus of transformer,
 * negative sign is opposite
 * type: integer
 * indexed
 */
#define TRANSFORMER_CONTROL "TRANSFORMER_CONTROL"

/**
 * Upper and lower limits of turns ratio or phase shift
 * type: float
 * indexed
 */
#define TRANSFORMER_RMA "TRANSFORMER_RMA"
#define TRANSFORMER_RMI "TRANSFORMER_RMI"

/**
 * Upper and lower limits of controlled volts, MW or MVAR
 * type: float
 * indexed
 */
#define TRANSFORMER_VMA "TRANSFORMER_VMA"
#define TRANSFORMER_VMI "TRANSFORMER_VMI"

/**
 * Number of available tap positions
 * type: integer
 * indexed
 */
#define TRANSFORMER_NTP "TRANSFORMER_NTP"

/**
 * Number of a transformer impedence table
 * type: integer
 * indexed
 */
#define TRANSFORMER_TAB "TRANSFORMER_TAB"

/**
 * Load drop compensation impedance for voltage controlling transformers
 * type: float
 * indexed
 */
#define TRANSFORMER_CR "TRANSFORMER_CR"
#define TRANSFORMER_CX "TRANSFORMER_CX"

/**
 * Winding connection angle in degrees
 * type: float
 * indexed
 */
#define TRANSFORMER_CNXA "TRANSFORMER_CNXA"

/**
 * Turns ratio increment
 * type: float
 * indexed
 */
#define TRANSFORMER_STEP "TRANSFORMER_STEP"

/**
 * Zero or number of a transformer impedence correction tabel (1-5)
 * type: integer
 * indexed
 */
#define TRANSFORMER_TABLE "TRANSFORMER_TABLE"

//  These transformer variables are not getting parsed for the V23 file format

/**
 * The winding one phase shift angle in degrees
 * Default value: 0.0
 * Type: real float
 * indexed
 */
#define TRANSFORMER_ANG1 "TRANSFORMER_ANG1"

/**
 * The winding data I/O code which defines the units in which
 * TRANSFORMER_WINDV1, and TRANSFORMER _WINDV2
 * are specified 
 * 1: off-nominal turns ratio in pu of winding bus base voltage
 * 2: winding voltage in kV.
 * Default value: 1
 * type: integer
 * indexed
 */
#define TRANSFORMER_CW "TRANSFORMER_CW"

/**
 * The impedance data I/O code defining the units in which R1-2, and X1-2 are specified
 * 1: for resistance and reactance in pu on system base quantities;
 * 2: for resistance and reactance in pu on a specified base MVA and winding bus base voltage
 * Default value: 1
 * type: integer
 * indexed
 */
#define TRANSFORMER_CZ "TRANSFORMER_CZ"

/**
 * The magnetizing admittance I/O code that defines the units in which TRANSFORMER_MAG1 and 
 * TRANSFORMER_MAG2 are specified
 * 1: for complex admittance in pu on system base quantities
 * 2: for no load loss in watts and exciting current in pu on winding one to two base MVA and
 *    nominal voltage
 * Default value: 1
 * type: integer
 * indexed
 */
#define TRANSFORMER_CM "TRANSFORMER_CM"

/**
 * The magnetizing conductance, in pu on system base quantities when TRANSFORMER_CM is 1; 
 * TRANSFORMER_MAG1 is the no load loss in watts when TRANSFORMER_CM is 2
 * type: real float
 * indexed
 */
#define TRANSFORMER_MAG1 "TRANSFORMER_MAG1"

/**
 * The magnetizing susceptance, in pu on system base quantities when CM is 1; 
 * TRANSFORMER_MAG2 is the exciting current in pu on winding one to two base MVA (SBASE1-2)
 * and nominal voltage (NOMV1) when TRANSFORMER_CM is 2
 * type: real float
 * indexed
 */
#define TRANSFORMER_MAG2 "TRANSFORMER_MAG2"

/**
 * The nonmetered end code of either 1 (for the winding one bus) or 2 (for the winding two bus).
 * Default value: 2
 * type: integer
 * indexed
 */
#define TRANSFORMER_NMETR "TRANSFORMER_NMETR"

/**
 * An alphanumeric identifier assigned to the transformer
 * type: string
 * indexed
 */
#define TRANSFORMER_NAME "TRANSFORMER_NAME"

/**
 * The initial transformer status, where 1 designates in-service and 0 designates out-of-service.
 * Default value: 1
 * type: integer
 * indexed
 */
#define TRANSFORMER_STATUS "TRANSFORMER_STATUS"

/**
 * Transformer owner number
 * type: integer
 * indexed
 */
#define TRANSFORMER_OWNER "TRANSFORMER_OWNER"

/**
 * The measured resistance of the transformer between the buses to which its first and second
 * windings are connected. 
 * Default value: 0.0
 * type: real float
 * indexed
 */
#define TRANSFORMER_R1_2 "TRANSFORMER_R1_2"

/**
 * The measured reactance of the transformer between the buses to which its
 * first and second
 * windings are connected.
 * type: real float
 * indexed
 */
#define TRANSFORMER_X1_2 "TRANSFORMER_X1_2"

/**
 * The winding one to two base MVA of the transformer
 * indexed
 * type: real float
 */
#define TRANSFORMER_SBASE1_2 "TRANSFORMER_SBASE1_2"

/**
 * Winding 1 and 2. Depends on CW
 * Default value: 1.0
 * type: real float
 * indexed
 */
#define TRANSFORMER_WINDV1 "TRANSFORMER_WINDV1"
#define TRANSFORMER_WINDV2 "TRANSFORMER_WINDV2"

/**
 * Nominal winding 1 and 2 voltage base.
 * Default value: 0.0
 * type: real float
 * indexed
 */
#define TRANSFORMER_NOMV1 "TRANSFORMER_NOMV1"
#define TRANSFORMER_NOMV2 "TRANSFORMER_NOMV2"


/**
 * Transformer control mode for winding 1 tap or phase shift
 * Default value: 0
 * type: integer
 * indexed
 */
#define TRANSFORMER_CODE1 "TRANSFORMER_CODE1"

/**
 * Bus number of the bus for which voltage is to be controlled by the
 * transformer turns ratio adjustment option
 * Default value: 0
 * type: integer
 * indexed
 */
#define TRANSFORMER_CONT1 "TRANSFORMER_CONT1"

/**
 * Number of tap positions available
 * Default value: 33
 * type: integer
 * indexed
 */
#define TRANSFORMER_NTP1 "TRANSFORMER_NTP1"

/**
 * The number of a transformer impedance correction table
 * Default value: 0
 * type: integer
 * indexed
 */
#define TRANSFORMER_TAB1 "TRANSFORMER_TAB1"

/**
 * Real and imaginary load drop compensation impedence for
 * voltage controlling transformers
 * Default value: 0.0
 * type: real float
 * indexed
 */
#define TRANSFORMER_CR1 "TRANSFORMER_CR1"
#define TRANSFORMER_CX1 "TRANSFORMER_CX1"

/**
 * Winding connection angle in degrees
 * Default value: 0.0
 * type: real float
 * indexed
 */
#define TRANSFORMER_CNXA1 "TRANSFORMER_CNXA1"

// Transformer Sequence Data
/**
 * Bus number to which another winding of transformer is connected
 * Default value: 0
 * type: integer
 * indexed
 */
#define TRANSFORMER_SEQ_K "TRANSFORMER_SEQ_K"

/**
 * Winding connection code
 * Default value: 4
 * type: integer
 * indexed
 */
#define TRANSFORMER_SEQ_CC "TRANSFORMER_SEQ_CC"

/**
 * Real part of zero sequence grounding impedance
 * Default value: 0.0
 * type: real float
 * indexed
 */
#define TRANSFORMER_SEQ_RG "TRANSFORMER_SEQ_RG"

/**
 * Imaginary part of zero sequence grounding impedance
 * Default value: 0.0
 * type: real float
 * indexed
 */
#define TRANSFORMER_SEQ_XG "TRANSFORMER_SEQ_XG"

/**
 * Real part of zero sequence impedance of two-winding transformer
 * type: real float
 * indexed
 */
#define TRANSFORMER_SEQ_R1 "TRANSFORMER_SEQ_R1"

/**
 * Imaginary part of zero sequence impedance of two-winding transformer
 * type: real float
 * indexed
 */
#define TRANSFORMER_SEQ_X1 "TRANSFORMER_SEQ_X1"

/**
 * Real zero sequence grounding impedance at the winding 2 side of an impedance
 * grounded two-winding transformer with connection code 8.
 * Default value: 0.0
 * type: real float
 * indexed
 */
#define TRANSFORMER_SEQ_RG2 "TRANSFORMER_SEQ_RG2"

/**
 * Imaginary zero sequence grounding impedance at the winding 2 side of an impedance
 * grounded two-winding transformer with connection code 8.
 * Default value: 0.0
 * type: real float
 * indexed
 */
#define TRANSFORMER_SEQ_XG2 "TRANSFORMER_SEQ_XG2"

/**
 * Real part winding two zero sequence impedance of a three-winding transformer
 * type: real float
 * indexed
 */
#define TRANSFORMER_SEQ_R2 "TRANSFORMER_SEQ_R2"

/**
 * Imaginary part winding two zero sequence impedance of a three-winding transformer
 * type: real float
 * indexed
 */
#define TRANSFORMER_SEQ_X2 "TRANSFORMER_SEQ_X2"

/**
 * Real part winding two threeo sequence impedance of a three-winding transformer
 * type: real float
 * indexed
 */
#define TRANSFORMER_SEQ_R3 "TRANSFORMER_SEQ_R3"

/**
 * Imaginary part winding three zero sequence impedance of a three-winding transformer
 * type: real float
 * indexed
 */
#define TRANSFORMER_SEQ_X3 "TRANSFORMER_SEQ_X3"

#endif /* _TRANSFORMER_VAR_HPP_ */
