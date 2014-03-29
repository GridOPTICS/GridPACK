/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all parameters that can be read in from
 * PTI format files. Each parameter has a corresponding macro that can be used
 * as a unique string to identify the parameter. The use of macros instead of
 * using strings directly will provide extra safety by forcing compiler errors
 * in the case of typos or spelling mistakes.
 */

/**
 *  Variables that can be associated more than once for a bus or a branch can be
 *  indexed by an integer to distinguish different instances. For example,
 *  multiple generators can be associated with a bus and the variables
 *  associated with each generator. The variables that have an associated index
 *  are denoted with the keyword "indexed".
 */

#ifndef DICTIONARY_HPP_
#define DICTIONARY_HPP_

// CASE DATA
/**
 * 0: base case
 * 1: add information to existing case
 * type: integer
 */
#define CASE_ID "CASE_ID"

/**
 * System base MVS. 
 * Default value is 100.0 MVA
 * type: real float
 */
#define CASE_SBASE "CASE_SBASE"

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
 * Active component of the shunt admittance to ground, entered in MW. 
 * Default value is 0.0
 * type: real float
 */
#define BUS_SHUNT_GL "BUS_SHUNT_GL"

/**
 * Reactive component of shunt admittance to ground, entered in Mvar. 
 * positive for a capacitor and negative for a reactor
 * Default value is 0.0
 * type: real float
 */
#define BUS_SHUNT_BL "BUS_SHUNT_BL"

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
 * Owner number
 * type: integer
 */
#define BUS_OWNER "BUS_OWNER"


// LOAD DATA

/**
 * The Bus number to which the load is connected
 * type: integer
 */
#define LOAD_BUSNUMBER "LOAD_BUSNUMBER"

/**
 * Non-blank alphanumeric identifier to distinguish different loads connected to the same bus.
 * Default value: ’1’
 * type: string
 */
#define LOAD_ID "LOAD_ID"

/**
 * Load status
 *  1: in-service
 *  0: out-of-service
 * Default value is 1
 * type: integer
 */
#define LOAD_STATUS "LOAD_STATUS"

/**
 * Area to which the load is assigned
 * type: integer
 */
#define LOAD_AREA "LOAD_AREA"

/**
 * Zone to which the load is assigned
 * type: integer
 */
#define LOAD_ZONE "LOAD_ZONE"

/**
 * Active power component of constant power load; entered in MW
 * type: real float
 */
#define LOAD_PL "LOAD_PL"

/**
 * Reactive power component of constant power load; entered in MVar
 * type: real float
 */
#define LOAD_QL "LOAD_QL"

/**
 * Active power component of constant current load; entered in MW 
 * type: real float
 */
#define LOAD_IP "LOAD_IP"

/**
 * Reactive power component of constant current load; entered in Mvar
 * type: real float
 */
#define LOAD_IQ "LOAD_IQ"

/**
 * Active power component of constant admittance load; entered in MW
 * type: real float
 */
#define LOAD_YP "LOAD_YP"

/**
 * Reactive power component of constant admittance load; entered in MVar
 * type: real float
 */
#define LOAD_YQ "LOAD_YQ"

/**
 * Owner to which the load is assigned
 * type: integer
 */
#define LOAD_OWNER "LOAD_OWNER"


// GENERATOR DATA
/**
 * Number of generators on a bus
 * type: integer
 */
#define GENERATOR_NUMBER "GENERATOR_NUMBER"

/**
 * Bus number to which the generator is connected
 * type: integer
 * indexed
 */
#define GENERATOR_BUSNUMBER "GENERATOR_BUSNUMBER"

/**
 * Non-blank alphanumeric machine identifier, used to distinguish among multiple machines connected to the same bus	
 * type: string
 * indexed
 */
#define GENERATOR_ID "GENERATOR_ID"

/**
 * Generator active power output, entered in MW	
 * type: real float
 * indexed
 */
#define GENERATOR_PG "GENERATOR_PG"

/**
 * Generator reactive power output, entered in MVar
 * type: real float
 * indexed
 */
#define GENERATOR_QG "GENERATOR_QG"

/**
 * Maximum generator reactive power output; entered in Mvar
 * type: real float
 * indexed
 */
#define GENERATOR_QMAX "GENERATOR_QMAX"

/**
 * Minimum generator reactive power output; entered in Mvar
 * type: real float
 * indexed
 */
#define GENERATOR_QMIN "GENERATOR_QMIN"

/**
 * Voltage setpoint; entered in pu
 * type: real float
 * indexed
 */
#define GENERATOR_VS "GENERATOR_VS"

/**
 * Bus number of a remote type 1 or 2 bus whose voltage is to be regulated by this plant to the
 * value specified by GENERATOR_VS
 * type: integer
 * indexed
 */
#define GENERATOR_IREG "GENERATOR_IREG"

/**
 * Total MVA base of the units represented by this machine; entered in MVA. 
 * type: real float
 * indexed
 */
#define GENERATOR_MBASE "GENERATOR_MBASE"

/**
 * Complex impedance, in pu on GENERATOR_MBASE base. 
 * type: complex
 * indexed
 */
#define GENERATOR_ZSOURCE "GENERATOR_ZSOURCE"

/**
 * Step-up transformer impedance; entered in pu on GENERATOR_MBASE base. 
 * type: complex
 * indexed
 */
#define GENERATOR_XTRAN "GENERATOR_XTRAN"

/**
 * Active and reactive components of Step-up transformer impedance, in pu on GENERATOR_MBASE base
 * type: real float
 * indexed
 */
#define GENERATOR_RT "GENERATOR_RT"
#define GENERATOR_XT "GENERATOR_XT"

/**
 * Step-up transformer off-nominal turns ratio; entered in pu
 * type: real float
 * indexed
 */
#define GENERATOR_GTAP "GENERATOR_GTAP"

/**
 * Initial machine status
 * 1: in-service
 * 0: out-of-service
 * type: integer
 * indexed
 */
#define GENERATOR_STAT "GENERATOR_STAT"

/**
 * Percent of the total Mvar required to hold the voltage at the bus controlled by bus 
 * that are to be contributed by the generation. It must be positive
 * type: real float
 * indexed
 */
#define GENERATOR_RMPCT "GENERATOR_RMPCT"

/**
 * Maximum generator active power output; entered in MW
 * type: real float
 * indexed
 */
#define GENERATOR_PMAX "GENERATOR_PMAX"

/**
 * Minimum generator active power output; entered in MW
 * type: real float
 * indexed
 */
#define GENERATOR_PMIN "GENERATOR_PMIN"

/**
 * Generator owner number	
 * type: integer
 * indexed
 */
#define GENERATOR_OWNER "GENERATOR_OWNER"

/**
 * Generator reactance
 * type: real float
 * indexed
 */
#define GENERATOR_REACTANCE "GENERATOR_REACTANCE"

/**
 * Generator resistance
 * type: real float
 * indexed
 */
#define GENERATOR_RESISTANCE "GENERATOR_RESISTANCE"

/**
 * Generator transient reactance
 * type: real float
 * indexed
 */
#define GENERATOR_TRANSIENT_REACTANCE "GENERATOR_TRANSIENT_REACTANCE"

/**
 * Generator subtransient reactance
 * type: real float
 * indexed
 */
#define GENERATOR_SUBTRANSIENT_REACTANCE "GENERATOR_SUBTRANSIENT_REACTANCE"

/**
 * Generator inertia constant
 * type: real float
 * indexed
 */
#define GENERATOR_INERTIA_CONSTANT_H "GENERATOR_INERTIA_CONSTANT_H"

/**
 * Generator damping coefficient
 * type: real float
 * indexed
 */
#define GENERATOR_DAMPING_COEFFICIENT_0 "GENERATOR_DAMPING_COEFFICENT_0"

// BRANCH DATA
/**
 * Global index used to sort branches into a fixed order
 * type: integer
 */
#define BRANCH_INDEX "BRANCH_INDEX"
/**
 * Branch “from bus”	
 * type: integer
 */
#define BRANCH_FROMBUS "BRANCH_FROMBUS"

/**
 * Branch “to bus”	
 * type: integer
 */
#define BRANCH_TOBUS "BRANCH_TOBUS"

/**
 * Number of transmission elements on branch
 * type: integer
 */
#define BRANCH_NUM_ELEMENTS "BRANCH_NUM_ELEMENTS"

/**
 * Non-blank alphanumeric branch circuit identifier
 * type: string
 * indexed
 */
#define BRANCH_CKT "BRANCH_CKT"

/**
 * Branch resistance; entered in pu
 * type: real float
 * indexed
 */
#define BRANCH_R "BRANCH_R"

/**
 * Branch reactance; entered in pu. 
 * type: real float
 * indexed
 */
#define BRANCH_X "BRANCH_X"

/**
 * Total branch charging susceptance; entered in pu
 * type: real float
 * indexed
 */
#define BRANCH_B "BRANCH_B"

/**
 * First current rating; entered in MVA
 * type: real float
 * indexed
 */
#define BRANCH_RATING_A "BRANCH_RATING_A"

/**
 * Second current rating; entered in MVA
 * type: real float
 * indexed
 */
#define BRANCH_RATING_B "BRANCH_RATING_B"

/**
 * Third current rating; entered in MVA
 * type: real float
 * indexed
 */
#define BRANCH_RATING_C "BRANCH_RATING_C"

/**
 * Transformer tap ratio in PTI 23 version
 * type: real float
 * indexed
 */
#define BRANCH_TAP "BRANCH_TAP"

/**
 * Transformer shift in PTI 23 version
 * type: real float
 * indexed
 */
#define BRANCH_SHIFT "BRANCH_SHIFT"

/**
 * Real part of admittance of the line shunt at the “from bus” end of the branch
 * type: real float
 * indexed
 */
#define BRANCH_SHUNT_ADMTTNC_G1 "BRANCH_SHUNT_ADMTTNC_G1"

/**
 * Imaginary part of admittance of the line shunt at the “from bus” end of the branch
 * type: real float
 * indexed
 */
#define BRANCH_SHUNT_ADMTTNC_B1 "BRANCH_SHUNT_ADMTTNC_B1"

/**
 * Real part of admittance of the line shunt at the “to bus” end of the branch
 * type: real float
 * indexed
 */
#define BRANCH_SHUNT_ADMTTNC_G2 "BRANCH_SHUNT_ADMTTNC_G2"

/**
 * Imaginary part of admittance of the line shunt at the “to bus” end of the branch
 * type: real float
 * indexed
 */
#define BRANCH_SHUNT_ADMTTNC_B2 "BRANCH_SHUNT_ADMTTNC_B2"


/**
 * Initial branch status
 * 1: in-service
 * 0: out-of-service
 * type: integer
 * indexed
 */
#define BRANCH_STATUS "BRANCH_STATUS"


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
 * The winding data I/O code which defines the units in which TRANSFORMER_WINDV1, and TRANSFORMER _WINDV2
 *  are specified 
 * 1: off-nominal turns ratio in pu of winding bus base voltage
 * 2: winding voltage in kV.
 * Default value: 1
 * type: integer
 */
#define TRANSFORMER_CW "TRANSFORMER_CW"

/**
 * The impedance data I/O code defining the units in which R1-2, and X1-2 are specified
 * 1: for resistance and reactance in pu on system base quantities;
 * 2: for resistance and reactance in pu on a specified base MVA and winding bus base voltage
 * Default value: 1
 * type: integer
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
 */
#define TRANSFORMER_CM "TRANSFORMER_CM"

/**
 * The magnetizing conductance, in pu on system base quantities when TRANSFORMER_CM is 1; 
 * TRANSFORMER_MAG1 is the no load loss in watts when TRANSFORMER_CM is 2
 * type: real float
 */
#define TRANSFORMER_MAG1 "TRANSFORMER_MAG1"

/**
 * The magnetizing susceptance, in pu on system base quantities when CM is 1; 
 * TRANSFORMER_MAG2 is the exciting current in pu on winding one to two base MVA (SBASE1-2)
 * and nominal voltage (NOMV1) when TRANSFORMER_CM is 2
 * type: real float
 */
#define TRANSFORMER_MAG2 "TRANSFORMER_MAG2"

/**
 * The nonmetered end code of either 1 (for the winding one bus) or 2 (for the winding two bus).
 * Default value: 2
 * type: integer
 */
#define TRANSFORMER_NMETR "TRANSFORMER_NMETR"

/**
 * An alphanumeric identifier assigned to the transformer
 * type: string
 */
#define TRANSFORMER_NAME "TRANSFORMER_NAME"

/**
 * The initial transformer status, where 1 designates in-service and 0 designates out-of-service.
 * Default value: 1
 * type: integer
 */
#define TRANSFORMER_STATUS "TRANSFORMER_STATUS"

/**
 * Transformer owner number
 * type: integer
 */
#define TRANSFORMER_OWNER "TRANSFORMER_OWNER"

/**
 * The measured resistance of the transformer between the buses to which its first and second
 * windings are connected. 
 * Default value: 0.0
 * type: real float
 */
#define TRANSFORMER_R1_2 "TRANSFORMER_R1_2"

/**
 * The measured reactance of the transformer between the buses to which its first and second
 * windings are connected.
 * type: real float
 */
#define TRANSFORMER_X1_2 "TRANSFORMER_X1_2"

/**
 * The winding one to two base MVA of the transformer
 * type: real float
 */
#define TRANSFORMER_SBASE1_2 "TRANSFORMER_SBASE1_2"


// AREA DATA
/**
 * Area number
 * type: integer
 */
#define AREAINTG_NUMBER "AREAINTG_NUMBER"

/**
 * Bus number of the area slack bus
 * type: integer
 */
#define AREAINTG_ISW "AREAINTG_ISW"

/**
 * Desired net interchange leaving the area (export); entered in MW
 * type: real float
 */
#define AREAINTG_PDES "AREAINTG_PDES"

/**
 * Interchange tolerance bandwidth; entered in MW
 * type: real float
 */
#define AREAINTG_PTOL "AREAINTG_PTOL"

/**
 * Area name
 * type: string
 */
#define AREAINTG_NAME "AREAINTG_NAME"


// SWITCHED DATA

/**
 * Bus number to which the shunt is connected
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
 * When SHUNT_MODSW is 1 or 2, the controlled voltage upper limit; entered in pu.
 * When SHUNT_MODSW is 3, 4 or 5, the controlled reactive power range upper limit;
 * entered in pu of the total reactive power range of the controlled voltage controlling device.
 * SHUNT_VSWHI = 1.0 by default.
 * type: real float
 */
#define SHUNT_VSWHI "SHUNT_VSWHI"

/**
 * When SHUNT_MODSW is 1 or 2, the controlled voltage lower limit; entered in pu.
 * When SHUNT_MODSW is 3, 4 or 5, the controlled reactive power range lower limit;
 * entered in pu of the total reactive power range of the controlled voltage controlling
 * device. SHUNT_VSWLO = 1.0 by default
 * type: real float
 */
#define SHUNT_VSWLO "SHUNT_VSWLO"

/**
 * Bus number of the bus whose voltage or connected equipment reactive power output
 * is controlled by this switched shunt
 * type: integer
 */
#define SHUNT_SWREM "SHUNT_SWREM"

/**
 * Initial switched shunt admittance; entered in Mvar at unity voltage
 * Default value: 0.0
 * type: real float
 */
#define SHUNT_BINIT "SHUNT_BINIT"

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


// TRANSFORMER IMPEDANCE CORRECTION
/**
 * Impedance correction table number
 * type: integer
 */
#define XFMR_CORR_TABLE_NUMBER "XFMR_CORR_TABLE_NUMBER"

/**
 * Either off-nominal turns ratio in pu or phase shift angle in degrees
 * type: real float
 */
#define XFMR_CORR_TABLE_Ti "XFMR_CORR_TABLE_Ti"

/**
 * Scaling factor by which transformer nominal impedance is to be multiplied to obtain
 * the actual transformer impedance for the corresponding "Ti"
 * type: real float
 */
#define XFMR_CORR_TABLE_Fi "XFMR_CORR_TABLE_Fi"


// MULTISECTION LINE GROUPING 
/**
 * “From bus” number of multisection line
 * type: integer
 */
#define MULTI_SEC_LINE_FROMBUS "MULTI_SEC_LINE_FROMBUS"

/**
 * “To bus” number. It is entered as a negative number or with a minus sign before the
 * first character of the extended bus name to designate it as the metered end; otherwise,
 * MULTI_SEC_LINE_FROMBUS is assumed to be the metered end
 * type: integer
 */
#define MULTI_SEC_LINE_TOBUS "MULTI_SEC_LINE_TOBUS"

/**
 * Two-character upper case alphanumeric multisection line grouping identifier.
 * type: string
 */
#define MULTI_SEC_LINE_ID "MULTI_SEC_LINE_ID"

/**
 * Bus numbers of the "dummy buses" connected by the branches that comprise this multisection
 * line grouping. No defaults allowed.
 * type: integer
 */
#define MULTI_SEC_LINE_DUMi "MULTI_SEC_LINE_DUMi"


// ZONE DATA
/**
 * Zone Number
 * type: integer
 */
#define ZONE_NUMBER "ZONE_NUMBER"

/**
 * Zone Name
 * type: string
 */
#define ZONE_NAME "ZONE_NAME"


// INTERAREA TRANSFER
/**
 * From area number of interarea transfer
 * type: integer
 */
#define INTERAREA_TRANSFER_FROM "INTERAREA_TRANSFER_FROM"

/**
 * To area number of interarea transfer
 * type: integer
 */
#define INTERAREA_TRANSFER_TO "INTERAREA_TRANSFER_TO"

/**
 * Single-character (0 through 9 or A through Z) upper case interarea transfer identifier
 * used to distinguish among multiple transfers
 * type: character
 */
#define INTERAREA_TRANSFER_TRID "INTERAREA_TRANSFER_TRID"

/**
 * MW comprising this transfer
 * type: real float
 */
#define INTERAREA_TRANSFER_PTRAN "INTERAREA_TRANSFER_PTRAN"


// OWNER
/**
 * Owner number
 * type: integer
 */
#define OWNER_NUMBER "OWNER_NUMBER"

/**
 * Owner name
 * type: integer
 */
#define OWNER_NAME "OWNER_NAME"



#endif /* DICTIONARY_HPP_ */
