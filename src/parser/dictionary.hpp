/*
 * This file contains definitions for all parameters that can be read in from
 * PTI format files. Each parameter has a corresponding macro that can be used
 * as a unique string to identify the parameter. The use of macros instead of
 * using strings directly will provide extra safety by forcing compiler errors
 * in the case of typos or spelling mistakes.
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
 * Default value is 100.0
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
 * Default value is 0.0
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
 * entered as part of load data.
 * Default value is 0.0
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


// LOAD DATA

/**
 * Bus number to which the load is connected
 * type: integer
 */
#define LOAD_BUSNUMBER "LOAD_BUSNUMBER"

/**
 * One- or two-character uppercase nonblank alphanumeric load identifier used to distinguish 
 * among multiple loads connected to the same bus.
 * Default value is ’1’
 * type: string
 */
#define LOAD_ID "LOAD_ID"

/**
 * Initial load status
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
 * Active power component of constant MVA load; entered in MW
 * type: real float
 */
#define LOAD_PL "LOAD_PL"

/**
 * Reactive power component of constant MVA load; entered in MVar
 * type: real float
 */
#define LOAD_QL "LOAD_QL"

/**
 * Active power component of constant current load; entered in MW at one per unit voltage
 * type: real float
 */
#define LOAD_IP "LOAD_IP"

/**
 * Reactive power component of constant current load; entered in Mvar at one per unit voltage
 * type: real float
 */
#define LOAD_IQ "LOAD_IQ"

/**
 * Active power component of constant admittance load; entered in MW at one per unit voltage
 * type: real float
 */
#define LOAD_YP "LOAD_YP"

/**
 * Reactive power component of constant admittance load; entered in MVar at one per unit voltage
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
 * Bus number to which the generator is connected
 * type: integer
 */
#define GENERATOR_BUSNUMBER "GENERATOR_BUSNUMBER"

/**
 * One- or two-character uppercase nonblank alphanumeric machine identifier 
 * used to distinguish among multiple machines connected to the same bus	
 * type: string
 */
#define GENERATOR_ID "GENERATOR_ID"

/**
 * Generator active power output, entered in MW	
 * type: real float
 */
#define GENERATOR_PG "GENERATOR_PG"

/**
 * Generator reactive power output, entered in MVar
 * type: real float
 */
#define GENERATOR_QG "GENERATOR_QG"

/**
 * Maximum generator reactive power output; entered in Mvar
 * type: real float
 */
#define GENERATOR_QMAX "GENERATOR_QMAX"

/**
 * Minimum generator reactive power output; entered in Mvar
 * type: real float
 */
#define GENERATOR_QMIN "GENERATOR_QMIN"

/**
 * Regulated voltage setpoint; entered in pu
 * type: real float
 */
#define GENERATOR_VS "GENERATOR_VS"

/**
 * Bus number of a remote type 1 or 2 bus whose voltage is to be regulated by this plant to the
 * value specified by GENERATOR_VS
 * type: integer
 */
#define GENERATOR_IREG "GENERATOR_IREG"

/**
 * Total MVA base of the units represented by this machine; entered in MVA. 
 * This quantity is not needed in normal power flow and equivalent construction work,
 * but is required for switching studies, fault analysis, and dynamic simulation.
 * type: real float
 */
#define GENERATOR_MBASE "GENERATOR_MBASE"

/**
 * Complex machine impedance entered in pu on GENERATOR_MBASE base. 
 * This data is not needed in normal power flow and equivalent construction work
 * but is required for switching studies, fault analysis, and dynamic simulation. 
 * For dynamic simulation, this impedance must be set equal to the subtransient impedance 
 * for those generators to be modeled by subtransient level machine models, 
 * and to transient impedance for those to be modeled by classical or transient level models
 * type: complex
 */
#define GENERATOR_ZSORCE "GENERATOR_ZSORCE"

/**
 * Step-up transformer impedance; entered in pu on GENERATOR_MBASE base. 
 * It should be entered as zero if the step-up transformer is explicitly modeled
 * type: complex
 */
#define GENERATOR_XTRAN "GENERATOR_XTRAN"

/**
 * Step-up transformer off-nominal turns ratio; entered in pu
 * type: real float
 */
#define GENERATOR_GTAP "GENERATOR_GTAP"

/**
 * Initial machine status
 * 1: in-service
 * 0: out-of-service
 * type: integer
 */
#define GENERATOR_STAT "GENERATOR_STAT"

/**
 * Percent of the total Mvar required to hold the voltage at the bus controlled by bus 
 * that are to be contributed by the generation. It must be positive
 * type: real float
 */
#define GENERATOR_RMPCT "GENERATOR_RMPCT"

/**
 * Maximum generator active power output; entered in MW
 * type: real float
 */
#define GENERATOR_PMAX "GENERATOR_PMAX"

/**
 * Minimum generator active power output; entered in MW
 * type: real float
 */
#define GENERATOR_PMIN "GENERATOR_PMIN"

/**
 * Generator owner number	
 * type: integer
 */
#define GENERATOR_OWNER "GENERATOR_OWNER"


// BRANCH DATA
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
 * One- or two-character uppercase nonblank alphanumeric branch circuit identifier
 * type: string
 */
#define BRANCH_CKT "BRANCH_CKT"

/**
 * Branch resistance; entered in pu
 * type: real float
 */
#define BRANCH_R "BRANCH_R"

/**
 * Branch reactance; entered in pu. A nonzero value of X must be entered for each branch
 * type: real float
 */
#define BRANCH_X "BRANCH_X"

/**
 * Total branch charging susceptance; entered in pu
 * type: real float
 */
#define BRANCH_B "BRANCH_B"

/**
 * First current rating; entered in MVA
 * type: real float
 */
#define BRANCH_RATING_A "BRANCH_RATING_A"

/**
 * Second current rating; entered in MVA
 * type: real float
 */
#define BRANCH_RATING_B "BRANCH_RATING_B"

/**
 * Third current rating; entered in MVA
 * type: real float
 */
#define BRANCH_RATING_C "BRANCH_RATING_C"

/**
 * Real part of admittance of the line shunt at the “from bus” end of the branch
 * type: real float
 */
#define BRANCH_SHUNT_ADMTTNC_G1 "BRANCH_SHUNT_ADMTTNC_G1"

/**
 * Imaginary part of admittance of the line shunt at the “from bus” end of the branch
 * type: real float
 */
#define BRANCH_SHUNT_ADMTTNC_B1 "BRANCH_SHUNT_ADMTTNC_B1"

/**
 * Real part of admittance of the line shunt at the “to bus” end of the branch
 * type: real float
 */
#define BRANCH_SHUNT_ADMTTNC_G2 "BRANCH_SHUNT_ADMTTNC_G2"

/**
 * Imaginary part of admittance of the line shunt at the “to bus” end of the branch
 * type: real float
 */
#define BRANCH_SHUNT_ADMTTNC_B2 "BRANCH_SHUNT_ADMTTNC_B2"

/**
 * Initial branch status
 * 1: in-service
 * 0: out-of-service
 * type: integer
 */
#define BRANCH_STATUS "BRANCH_STATUS"

/**
 * Line length; entered in user-selected units
 * type: real float
 */
#define BRANCH_LENGTH "BRANCH_LENGTH"

/**
 * Branch owner number
 * type: integer
 */
#define BRANCH_OWNER "BRANCH_OWNER"


// TRANSFORMER DATA 
/**
 * Bus number to which the first winding is connected
 * type: integer
 */
#define TRANSFORMER_BUS1 "TRANSFORMER_BUS1"

/**
 * Bus number to which the second winding is connected
 * type: integer
 */
#define TRANSFORMER_BUS2 "TRANSFORMER_BUS2"

/**
 * Bus number to which the third winding is connected. Zero is used for 2-winding transformers
 * type: integer
 */
#define TRANSFORMER_BUS3 "TRANSFORMER_BUS3"

/**
 * One- or two-character uppercase nonblank alphanumeric transformer circuit identifier
 * type: string
 */
#define TRANSFORMER_CKT "TRANSFORMER_CKT"

/**
 * The winding data I/O code which defines the units in which TRANSFORMER_WINDV1, TRANSFORMER _WINDV2
 * and TRANSFORMER _WINDV3 are specified (the units of RMAi and RMIi are also governed by 
 * TRANSFORMER_CW when |TRANSFORMER_CODi| is 1 or 2)
 * 1: off-nominal turns ratio in pu of winding bus base voltage
 * 2: winding voltage in kV.
 * Default value: 1
 * type: integer
 */
#define TRANSFORMER_CW "TRANSFORMER_CW"

/**
 * The impedance data I/O code that defines the units in which R1-2, X1-2, R2-3, X2-3, 
 * R3-1 and X3-1 are specified
 * 1: for resistance and reactance in pu on system base quantities;
 * 2: for resistance and reactance in pu on a specified base MVA and winding bus base voltage
 * 3: for transformer load loss in watts and impedance magnitude in pu on a specified base MVA
 *    and winding bus base voltage.
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
 * In addition, for a three-winding transformer, 3 (for the winding three bus) is a valid
 * specification
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
 * In addition, for a three-winding transformer, 2 designates that only winding two is out-of-service,
 * 3 indicates that only winding three is out-of-service, and 4 indicates that only winding
 * one is out-of-service, with the remaining windings in-service.
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
 * windings are connected. When TRANSFORMER_CZ is 1, they are the resistance and reactance,
 * respectively, in pu on system base quantities; when TRANSFORMER_CZ is 2, they are the resistance
 * and reactance, respectively, in pu on winding one to two base MVA (TRANSFORMER_SBASE1-2)
 * and winding one bus base voltage; when TRANSFORMER_CZ is 3, TRANSFORMER_R1-2 is the load
 * loss in watts, and TRANSFORMER_X1-2 is the impedance magnitude in pu on winding one to two
 * base MVA (SBASE1-2) and winding one bus base voltage.
 * Default value: 0.0
 * type: real float
 */
#define TRANSFORMER_R1-2 "TRANSFORMER_R1-2"

/**
 * The measured reactance of the transformer between the buses to which its first and second
 * windings are connected.
 * no default is allowed for TRANSFORMER_X1-2.
 * type: real float
 */
#define TRANSFORMER_X1-2 "TRANSFORMER_X1-2"

/**
 * The winding one to two base MVA of the transformer
 * type: real float
 */
#define TRANSFORMER_SBASE1-2 "TRANSFORMER_SBASE1-2"

/**
 * The measured resistance of a three-winding transformer between the buses to which its second
 * and third windings are connected; ignored for a two-winding transformer
 * type: real float
 */
#define TRANSFORMER_R2-3 "TRANSFORMER_R2-3"

/**
 * The measured reactance of a three-winding transformer between the buses to which its second
 * and third windings are connected; ignored for a two-winding transformer
 * type: real float
 */
#define TRANSFORMER_X2-3 "TRANSFORMER_X2-3"

/**
 * The winding two to three base MVA of a three-winding transformer; ignored for a two-winding
 * transformer
 * type: real float
 */
#define TRANSFORMER_SBASE2-3 "TRANSFORMER_SBASE2-3"

/**
 * The measured resistance of a three-winding transformer between the buses to which its third
 * and first windings are connected; ignored for a two-winding transformer
 * type: real float
 */
#define TRANSFORMER_R3-1 "TRANSFORMER_R3-1"

/**
 * The measured reactance of a three-winding transformer between the buses to which its third
 * and first windings are connected; ignored for a two-winding transformer
 * type: real float
 */
#define TRANSFORMER_X3-1 "TRANSFORMER_X3-1"

/**
 * The winding three to one base MVA of a three-winding transformer; ignored for a two-winding
 * transformer
 * type: real float
 */
#define TRANSFORMER_SBASE3-1 "TRANSFORMER_SBASE3-1"

/**
 * The voltage magnitude at the hidden "star point" bus; entered in pu
 * type: real float
 */
#define TRANSFORMER_VMSTAR "TRANSFORMER_VMSTAR"

/**
 * The bus voltage phase angle at the hidden "star point" bus; entered in degrees
 * type: real float
 */
#define TRANSFORMER_ANSTAR "TRANSFORMER_ANSTAR"

/**
 * The winding one off-nominal turns ratio in pu of winding one bus base voltage when
 * TRANSFORMER_CW is 1; TRANSFORMER_WINDV1 = 1.0 by default. TRANSFORMER_WINDV1 is the
 * actual winding one voltage in kV when TRANSFORMER_CW is 2;
 * type: real float
 */
#define TRANSFORMER_WINDV1 "TRANSFORMER_WINDV1"

/**
 * The nominal (rated) winding one voltage in kV, or zero to indicate that nominal winding
 * one voltage is to be taken as the base voltage of bus TRANSFORMER_BUS1. TRANSFORMER_NOMV1 is used only
 * in converting magnetizing data between per unit admittance values and physical units when
 * TRANSFORMER_CM is 2
 * Default value: 0.0
 * type: real float
 */
#define TRANSFORMER_NOMV1 "TRANSFORMER_NOMV1"

/**
 * The winding one phase shift angle in degrees. TRANSFORMER_ANG1 is positive for a positive
 * phase shift from the winding one side to the winding two side (for a two-winding transformer),
 * or from the winding one side to the "T" (or "star") point bus (for a three-winding transformer).
 * TRANSFORMER_ANG1 must be greater than -180.0 and less than or equal to +180.0.
 * Default value: 0.0
 * type: real float
 */
#define TRANSFORMER_ANG1 "TRANSFORMER_ANG1"

/**
 * The rating of the first winding, entered in MVA
 * type: real float
 */
#define TRANSFORMER_RATA1 "TRANSFORMER_RATA1"

/**
 * The rating of the second winding, entered in MVA
 * type: real float
 */
#define TRANSFORMER_RATB1 "TRANSFORMER_RATB1"

/**
 * The rating of the third winding, entered in MVA
 * type: real float
 */
#define TRANSFORMER_RATC1 "TRANSFORMER_RATC1"

/**
 * The transformer control mode for automatic adjustments of the winding one tap or phase
 * shift angle during power flow solutions
 *  0: for no control (fixed tap and phase shift)
 * ±1: for voltage control
 * ±2: for reactive power flow control
 * ±3: for active power flow control
 * ±4: for control of a dc line quantity (+4 is valid only for two-winding transformers).
 * If the control mode is entered as a positive number, automatic adjustment of this
 * transformer winding is enabled when the corresponding adjustment is activated during
 * power flow solutions; a negative control mode suppresses the automatic adjustment of
 * this transformer winding
 * Default value: 0
 * type: integer
 */
#define TRANSFORMER_COD1 "TRANSFORMER_COD1"

/**
 * The bus number, of the bus whose voltage is to be controlled by the transformer
 * turns ratio adjustment option of the power flow solution activities when TRANSFORMER_COD1 is 1.
 * TRANSFORMER_CONT1 should be nonzero only for voltage controlling transformer windings.
 * type: integer
 */
#define TRANSFORMER_CONT1 "TRANSFORMER_CONT1"

/**
 * The upper limit of either:
 * • Off-nominal turns ratio in pu of winding one bus base voltage when |TRANSFORMER_COD1|
 *   is 1 or 2 and TRANSFORMER_CW is 1; Default value: 1.1
 * • Actual winding one voltage in kV when |TRANSFORMER_COD1| is 1 or 2 and CW is 2.
 *   No default is allowed.
 * • Phase shift angle in degrees when |TRANSFORMER_COD1| is 3. No default is allowed. 
 * • Not used when |COD1| is 0 or 4; Default value: 1.1
 * type: real float
 */
#define TRANSFORMER_RMA1 "TRANSFORMER_RMA1"

/**
 * The lower limit of either:
 * • Off-nominal turns ratio in pu of winding one bus base voltage when |TRANSFORMER_COD1|
 *   is 1 or 2 and CW is 1; Default value: 0.9
 * • Actual winding one voltage in kV when |TRANSFORMER_COD1| is 1 or 2 and CW is 2.
 *   No default is allowed.
 * • Phase shift angle in degrees when |TRANSFORMER_COD1| is 3. No default is allowed. 
 * • Not used when |COD1| is 0 or 4; Default value: 0.9
 * type: real float
 */
#define TRANSFORMER_RMI1 "TRANSFORMER_RMI1"

/**
 * The upper limit of either:
 * • Voltage at the controlled bus (bus |TRANSFORMER_CONT1|) in pu when |TRANSFORMER_COD1| is 1.
 *   Default value: 1.1
 * • Reactive power flow into the transformer at the winding one bus end in Mvar when
 *   |TRANSFORMER_COD1| is 2. No default is allowed.
 * • Active power flow into the transformer at the winding one bus end in MW when
 *   |TRANSFORMER_COD1| is 3. No default is allowed.
 * • Not used when |COD1| is 0 or 4; Default value: 1.1
 * type: real float
 */
#define TRANSFORMER_VMA1 "TRANSFORMER_VMA1"

/**
 * The lower limit of either:
 *• Voltage at the controlled bus (bus |TRANSFORMER_CONT1|) in pu when |TRANSFORMER_COD1| is 1.
 *  Default value: 0.9
 *• Reactive power flow into the transformer at the winding one bus end in Mvar when
 *  |TRANSFORMER_COD1| is 2. No default is allowed.
 *• Active power flow into the transformer at the winding one bus end in MW when
 *  |TRANSFORMER_COD1| is 3. No default is allowed.
 *• Not used when |COD1| is 0 or 4; Default value: 0.9
 * type: real float
 */
#define TRANSFORMER_VMI1 "TRANSFORMER_VMI1"

/**
 * The number of tap positions available; used when TRANSFORMER_COD1 is 1 or 2. TRANSFORMER_NTP1
 * must be between 2 and 9999.
 * Default value: 33
 * type: integer
 */
#define TRANSFORMER_NTP1 "TRANSFORMER_NTP1"

/**
 * The number of a transformer impedance correction table if this transformer winding’s
 * impedance is to be a function of either off-nominal turns ratio or phase shift angle.
 * Default value: 0
 * type: integer
 */
#define TRANSFORMER_TAB1 "TRANSFORMER_TAB1"

/**
 * The load drop compensation resistance for voltage controlling transformers entered
 * in pu on system base quantities; used when TRANSFORMER_COD1 is 1.
 * type: real float
 */
#define TRANSFORMER_CR1 "TRANSFORMER_CR1"

/**
 * The load drop compensation reactance for voltage controlling transformers entered in pu
 * on system base quantities; used when TRANSFORMER_COD1 is 1.
 * type: real float
 */
#define TRANSFORMER_CX1 "TRANSFORMER_CX1"

/**
 * The winding two off-nominal turns ratio in pu of winding two bus base voltage when
 * TRANSFORMER_CW is 1; WINDV2 = 1.0 by default. TRANSFORMER_WINDV2 is the actual winding
 * two voltage in kV when TRANSFORMER_CW is 2; TRANSFORMER_WINDV2 is equal to the base voltage
 * of bus TRANSFORMER_BUS2 by default.
 * type: real float
 */
#define TRANSFORMER_WINDV2 "TRANSFORMER_WINDV2"

/**
 * The nominal (rated) winding two voltage in kV, or zero to indicate that nominal winding
 * two voltage is to be taken as the base voltage of bus TRANSFORMER_BUS2. TRANSFORMER_NOMV2
 * is present for information purposes only; it is not used in any of the calculations
 * for modeling the transformer.
 * Default value: 0.0
 * type: real float
 */
#define TRANSFORMER_NOMV2 "TRANSFORMER_NOMV2"

/**
 * The winding two phase shift angle in degrees; ignored for a two-winding transformer.
 * Positive for a positive phase shift from the winding two side to the "T" (or "star")
 * point bus. It must be greater than -180.0 and less than or equal to +180.0.
 * Default value: 0.0
 * type: real float
 */
#define TRANSFORMER_ANG2 "TRANSFORMER_ANG2"

/**
 * The second winding’s first ratings entered in MVA; ignored for a two-winding transformer.
 * Default value: 0.0
 * type: real float
 */
#define TRANSFORMER_RATA2 "TRANSFORMER_RATA2"

/**
 * The second winding’s second ratings entered in MVA; ignored for a two-winding transformer.
 * Default value: 0.0
 * type: real float
 */
#define TRANSFORMER_RATB2 "TRANSFORMER_RATB2"

/**
 * The second winding’s third ratings entered in MVA; ignored for a two-winding transformer.
 * Default value: 0.0
 * type: real float
 */
#define TRANSFORMER_RATC2 "TRANSFORMER_RATC2"

/**
 * The transformer control mode for automatic adjustments of the winding two tap or phase
 * shift angle during power flow solutions: 
 * 0: for no control (fixed tap and phase shift)
 * ±1 for voltage control
 * ±2 for reactive power flow control
 * ±3 for active power flow control
 * If the control mode is entered as a positive number, automatic adjustment of this
 * transformer winding is enabled when the corresponding adjustment is activated during
 * power flow solutions; a negative control mode suppresses the automatic adjustment of
 * this transformer winding.
 * Default value: 0
 * type: integer
 */
#define TRANSFORMER_COD2 "TRANSFORMER_COD2"

/**
 * The bus number of the bus whose voltage is to be controlled by the transformer turns
 * ratio adjustment option of the power flow solution activities when TRANSFORMER_COD2 is 1.
 * It should be nonzero only for voltage controlling transformer windings.
 * Default value: 0
 * type: integer
 */
#define TRANSFORMER_CONT2 "TRANSFORMER_CONT2"

/**
 * The upper limit of either:
 * • Off-nominal turns ratio in pu of winding two bus base voltage when
 *   |TRANSFORMER_COD2| is 1 or 2 and TRANSFORMER_CW is 1; Default value 1.1
 * • Actual winding two voltage in kV when |TRANSFORMER_COD2| is 1 or 2 and
 *   TRANSFORMER_CW is 2. No default is allowed.
 * • Phase shift angle in degrees when |TRANSFORMER_COD2| is 3. No default is allowed. 
 * • Not used when |TRANSFORMER_COD2| is 0; Default value: 1.1
 * type: real float
 */
#define TRANSFORMER_RMA2 "TRANSFORMER_RMA2"

/**
 * The lower limit of either:
 * • Off-nominal turns ratio in pu of winding two bus base voltage when
 *   |TRANSFORMER_COD2| is 1 or 2 and TRANSFORMER_CW is 1; Default value 0.9
 * • Actual winding two voltage in kV when |TRANSFORMER_COD2| is 1 or 2 and
 *   TRANSFORMER_CW is 2. No default is allowed.
 * • Phase shift angle in degrees when | TRANSFORMER_COD2| is 3. No default is allowed. 
 * • Not used when |TRANSFORMER_COD2| is 0; Default value: 0.9
 * type: real float
 */
#define TRANSFORMER_RMI2 "TRANSFORMER_RMI2"

/**
 * The upper limit, of either:
 * • Voltage at the controlled bus (bus |TRANSFORMER_CONT2|) in pu when
 *   |TRANSFORMER_COD2| is 1. Default value: 1.1
 * • Reactive power flow into the transformer at the winding two bus end in Mvar when
 *   |TRANSFORMER_COD2| is 2. No default is allowed.
 * • Active power flow into the transformer at the winding two bus end in MW when
 *   |TRANSFORMER_COD2| is 3. No default is allowed.
 * • Not used when | TRANSFORMER_COD2| is 0; Default value: 1.1
 * type: real float
 */
#define TRANSFORMER_VMA2 "TRANSFORMER_VMA2"

/**
 * The lower limit, of either:
 * • Voltage at the controlled bus (bus |TRANSFORMER_CONT2|) in pu when
 *   |TRANSFORMER_COD2| is 1. Default value: 0.9
 * • Reactive power flow into the transformer at the winding two bus end in Mvar when
 *   |TRANSFORMER_COD2| is 2. No default is allowed.
 * • Active power flow into the transformer at the winding two bus end in MW when
 *   |TRANSFORMER_COD2| is 3. No default is allowed.
 * • Not used when |TRANSFORMER_COD2| is 0; Default value: 0.9
 * type: real float
 */
#define TRANSFORMER_VMI2 "TRANSFORMER_VMI2"

/**
 * The number of tap positions available; used when TRANSFORMER_COD2 is 1 or 2.
 * Must be between 2 and 9999.
 * Default value: 33
 * type: integer
 */
#define TRANSFORMER_NTP2 "TRANSFORMER_NTP2"

/**
 * The number of a transformer impedance correction table if this transformer winding’s
 * impedance is to be a function of either off-nominal turns ratio or phase shift angle,
 * or 0 if no transformer impedance correction is to be applied to this transformer winding.
 * Default value: 0
 * type: integer
 */
#define TRANSFORMER_TAB2 "TRANSFORMER_TAB2"

/**
 * The load drop compensation resistance for voltage controlling transformers entered in pu
 * on system base quantities; used when TRANSFORMER_COD2 is 1.
 * Default value: 0
 * type: real float
 */
#define TRANSFORMER_CR2 "TRANSFORMER_CR2"

/**
 * The load drop compensation reactance for voltage controlling transformers entered in pu
 * on system base quantities; used when TRANSFORMER_ COD2 is 1.
 * Default value: 0
 * type: real float
 */
#define TRANSFORMER_CX2 "TRANSFORMER_CX2"

/**
 * The winding three off-nominal turns ratio in pu of winding three bus base voltage
 * when TRANSFORMER_CW is 1; Default value: 1.0. TRANSFORMER_WINDV3 is the actual
 * winding three voltage in kV when TRANSFORMER_CW is 2; Default value: TRANSFORMER_BUS3
 * type: real float
 */
#define TRANSFORMER_WINDV3 "TRANSFORMER_WINDV3"

/**
 * The nominal (rated) winding three voltage in kV, or zero to indicate that nominal
 * winding three voltage is to be taken as the base voltage of bus TRANSFORMER_BUS3.
 * It is present for information purposes only; it is not used in any of the calculations
 * for modeling the transformer.
 * Default value: 0.0
 * type: real float
 */
#define TRANSFORMER_NOMV3 "TRANSFORMER_NOMV3"

/**
 * The winding three phase shift angle in degrees. It is positive for a positive phase
 * shift from the winding three side to the "T" (or “star”) point bus. It must be greater
 * than -180.0 and less than or equal to +180.0
 * Default value: 0.0
 * type: real float
 */
#define TRANSFORMER_ANG3 "TRANSFORMER_ANG3"

/**
 * The third winding’s first rating entered in MVA
 * Default value: 0
 * type: real float
 */
#define TRANSFORMER_RATA3 "TRANSFORMER_RATA3"

/**
 * The third winding’s second rating entered in MVA
 * Default value: 0
 * type: real float
 */
#define TRANSFORMER_RATB3 "TRANSFORMER_RATB3"

/**
 * The third winding’s third rating entered in MVA
 * Default value: 0
 * type: real float
 */
#define TRANSFORMER_RATC3 "TRANSFORMER_RATC3"

/**
 * The transformer control mode for automatic adjustments of the winding three tap
 * or phase shift angle during power flow solutions:
 *  0: for no control (fixed tap and phase shift)
 * ±1: for voltage control
 * ±2: for reactive power flow control
 * ±3: for active power flow control.
 * If the control mode is entered as a positive number, automatic adjustment of
 * this transformer winding is enabled when the corresponding adjustment
 * is activated during power flow solutions; a negative control mode suppresses the
 * automatic adjustment of this transformer winding.
 * Default value: 0
 * type: integer
 */
#define TRANSFORMER_COD3 "TRANSFORMER_COD3"

/**
 * The bus number of the bus whose voltage is to be controlled by the transformer
 * turns ratio adjustment option of the power flow solution activities when
 * TRANSFORMER_COD3 is 1. It should be nonzero only for voltage controlling transformer
 * type: integer
 */
#define TRANSFORMER_CONT3 "TRANSFORMER_CONT3"

/**
 * The upper limit of either:
 *• Off-nominal turns ratio in pu of winding three bus base voltage when
 *  |TRANSFORMER_COD3| is 1 or 2 and TRANSFORMER_CW is 1; Default value 1.1
 *• Actual winding three voltage in kV when |TRANSFORMER_COD3| is 1 or 2 and
 *  TRANSFORMER_CW is 2. No default is allowed.
 *• Phase shift angle in degrees when |TRANSFORMER_COD3| is 3. No default is allowed. 
 *• Not used when |TRANSFORMER_COD3| is 0
 * Default value: 1.1
 * type: real float
 */
#define TRANSFORMER_RMA3 "TRANSFORMER_RMA3"

/**
 * The lower limit, of either:
 *• Off-nominal turns ratio in pu of winding three bus base voltage when
 *  |TRANSFORMER_COD3| is 1 or 2 and TRANSFORMER_CW is 1; Default value 0.9
 *• Actual winding three voltage in kV when |TRANSFORMER_COD3| is 1 or 2 and
 *  TRANSFORMER_CW is 2. No default is allowed.
 *• Phase shift angle in degrees when |TRANSFORMER_COD3| is 3. No default is allowed. 
 *• Not used when |TRANSFORMER_COD3| is 0; Default value: 0.9
 * type: real float
 */
#define TRANSFORMER_RMI3 "TRANSFORMER_RMI3"

/**
 * The upper limit of either:
 *• Voltage at the controlled bus (bus |TRANSFORMER_CONT3|) in pu when
 *  |TRANSFORMER_COD3| is 1. Default value: 1.1
 *• Reactive power flow into the transformer at the winding two bus end in Mvar when
 *  |TRANSFORMER_COD3| is 2. No default is allowed.
 *• Active power flow into the transformer at the winding two bus end in MW when
 *  |TRANSFORMER_COD3| is 3. No default is allowed.
 *• Not used when | TRANSFORMER_COD3| is 0; Default value: 1.1
 * type: real float
 */
#define TRANSFORMER_VMA3 "TRANSFORMER_VMA3"

/**
 * The lower limit of either:
 *• Voltage at the controlled bus (bus | TRANSFORMER_CONT3|) in pu when
 *  |TRANSFORMER_COD3| is 1. Default value: 0.9
 *• Reactive power flow into the transformer at the winding two bus end in Mvar when
 *  |TRANSFORMER_COD3| is 2. No default is allowed.
 *• Active power flow into the transformer at the winding two bus end in MW when
 *  |TRANSFORMER_COD3| is 3. No default is allowed.
 *• Not used when |TRANSFORMER_COD3| is 0; Default value: 0.9
 * type: real float
 */
#define TRANSFORMER_VMI3 "TRANSFORMER_VMI3"

/**
 * The number of tap positions available; used when TRANSFORMER_COD3 is 1 or 2.
 * It must be between 2 and 9999.
 * Default value: 33
 * type: integer
 */
#define TRANSFORMER_NTP3 "TRANSFORMER_NTP3"

/**
 * The number of a transformer impedance correction table if this transformer winding’s
 * impedance is to be a function of either off-nominal turns ratio or
 * phase shift angle, or 0 if no transformer impedance correction is to be applied
 * to this transformer winding.
 * Default value: 0
 * type: integer
 */
#define TRANSFORMER_TAB3 "TRANSFORMER_TAB3"

/**
 * The load drop compensation resistance for voltage controlling transformers entered
 * in pu on system base quantities; used when TRANSFORMER_COD3 is 1.
 * Default value 0
 * type: real float
 */
#define TRANSFORMER_CR3 "TRANSFORMER_CR3"

/**
 * The load drop compensation reactance for voltage controlling transformers entered
 * in pu on system base quantities; used when TRANSFORMER_COD3 is 1.
 * Default value 0
 * type: real float
 */
#define TRANSFORMER_CX3 "TRANSFORMER_CX3"


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
 * SHUNT_VSWHI is not used when SHUNT_MODSW is 0. SHUNT_VSWHI = 1.0 by default.
 * type: real float
 */
#define SHUNT_VSWHI "SHUNT_VSWHI"

/**
 * When SHUNT_MODSW is 1 or 2, the controlled voltage lower limit; entered in pu.
 * When SHUNT_MODSW is 3, 4 or 5, the controlled reactive power range lower limit;
 * entered in pu of the total reactive power range of the controlled voltage controlling
 * device. SHUNT_VSWLO is not used when SHUNT_MODSW is 0. SHUNT_VSWLO = 1.0 by default
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
 * Percent of the total Mvar required to hold the voltage at the bus controlled by
 * the bus that are to be contributed by this switched shunt; SHUNT_RMPCT must be positive
 * type: real float
 */
#define SHUNT_RMPCT "SHUNT_RMPCT"

/**
 * When SHUNT_MODSW is 4, the name of the VSC dc line whose converter bus is specified
 * in SHUNT_SWREM. SHUNT_RMIDNT is not used for other values of SHUNT_MODSW.
 * Default value: a blank name
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
 * Number of steps for block i. The first zero value of Ni or Bi is interpreted as
 * the end of the switched shunt blocks for bus I.
 * Default value: 0
 * type: integer
 */
#define SHUNT_Ni "SHUNT_Ni"

/**
 * Admittance increment for each of Ni steps in block i; entered in Mvar at unity voltage.
 * Default value: 0
 * type: real float
 */
#define SHUNT_Bi "SHUNT_Bi"


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
