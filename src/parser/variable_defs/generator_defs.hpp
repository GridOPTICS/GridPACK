/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all generator parameters that can be read
 * in from * PTI format files. Each parameter has a corresponding macro that can
 * be used as a unique string to identify the parameter. The use of macros
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


/**
 * Variables that have _CURRENT appended to them may be added to the data
 * collection during runtime and updated as the simulation
 * proceeds
 */

#ifndef _GENERATOR_VAR_HPP_
#define _GENERATOR_VAR_HPP_

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
 * Non-blank alphanumeric machine identifier, used to distinguish
 * among multiple machines connected to the same bus	
 * type: string
 * indexed
 */
#define GENERATOR_ID "GENERATOR_ID"

/**
 * Unique global index that runs sequentially over all generators
 * type: integer
 * indexed
 */
#define GENERATOR_GLOBAL_IDX "GENERATOR_GLOBAL_IDX"

/**
 * Generator active power output, entered in MW	
 * type: real float
 * indexed
 */
#define GENERATOR_PG "GENERATOR_PG"
#define GENERATOR_PG_CURRENT "GENERATOR_PG_CURRENT"

/**
 * Generator reactive power output, entered in MVar
 * type: real float
 * indexed
 */
#define GENERATOR_QG "GENERATOR_QG"
#define GENERATOR_QG_CURRENT "GENERATOR_QG_CURRENT"

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
#define GENERATOR_VS "GENERATOR_VS_CURRENT"

/**
 * Bus number of a remote type 1 or 2 bus whose voltage is to be regulated by this plant to the
 * value specified by GENERATOR_VS
 * type: integer
 * indexed
 */
#define GENERATOR_IREG "GENERATOR_IREG"

/**
 * Node number of bus IREG
 * type: integer
 * indexed
 */
#define GENERATOR_NREG "GENERATOR_NREG"

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
 * Base loaded flag
 *     0: normal, not base loaded
 *     1: down only, machine can only be scaled down
 *     2: neither up or down
 *     3: up only, machine can only be scaled up
 * type: integer
 * indexed
 */
#define GENERATOR_BASLOD "GENERATOR_BASLOD"

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

/**
 * Generator TDOP
 * type: real float
 * indexed
 */
#define GENERATOR_TDOP "GENERATOR_TDOP"

/**
 * Generator TDOPP
 * type: real float
 * indexed
 */
#define GENERATOR_TDOPP "GENERATOR_TDOPP"

/**
 * Generator TQOP
 * type: real float
 * indexed
 */
#define GENERATOR_TQOP "GENERATOR_TQOP"

/**
 * Generator TQOPP
 * type: real float
 * indexed
 */
#define GENERATOR_TQOPP "GENERATOR_TQOPP"

/**
 * Generator XD
 * type: real float
 * indexed
 */
#define GENERATOR_XD "GENERATOR_XD"

/**
 * Generator XQ
 * type: real float
 * indexed
 */
#define GENERATOR_XQ "GENERATOR_XQ"

/**
 * Generator XDP
 * type: real float
 * indexed
 */
#define GENERATOR_XDP "GENERATOR_XDP"

/**
 * Generator XQP
 * type: real float
 * indexed
 */
#define GENERATOR_XQP "GENERATOR_XQP"

/**
 * Generator XDPP
 * type: real float
 * indexed
 */
#define GENERATOR_XDPP "GENERATOR_XDPP"

/**
 * Generator XL
 * type: real float
 * indexed
 */
#define GENERATOR_XL "GENERATOR_XL"

/**
 * Generator S1
 * type: real float
 * indexed
 */
#define GENERATOR_S1 "GENERATOR_S1"

/**
 * Generator S12
 * type: real float
 * indexed
 */
#define GENERATOR_S12 "GENERATOR_S12"

/**
 * Generator owner 1
 * type: real float
 * indexed
 */
#define GENERATOR_OWNER1 "GENERATOR_OWNER1"

/**
 * Generator owner 2
 * type: integer
 * indexed
 */
#define GENERATOR_OWNER2 "GENERATOR_OWNER2"

/**
 * Generator owner 3
 * type: integer
 * indexed
 */
#define GENERATOR_OWNER3 "GENERATOR_OWNER3"

/**
 * Generator owner 4
 * type: integer
 * indexed
 */
#define GENERATOR_OWNER4 "GENERATOR_OWNER4"

/**
 * Generator owner 1 fraction
 * type: integer
 * indexed
 */
#define GENERATOR_OFRAC1 "GENERATOR_OFRAC1"

/**
 * Generator owner 2 fraction
 * type: real float
 * indexed
 */
#define GENERATOR_OFRAC2 "GENERATOR_OFRAC2"

/**
 * Generator owner 3 fraction
 * type: real float
 * indexed
 */
#define GENERATOR_OFRAC3 "GENERATOR_OFRAC3"

/**
 * Generator owner 4 fraction
 * type: real float
 * indexed
 */
#define GENERATOR_OFRAC4 "GENERATOR_OFRAC4"

/**
 * Generator wind mode (distinguish wind generators)
 * type: integer
 * indexed
 */
#define GENERATOR_WMOD "GENERATOR_WMOD"

/**
 * Wind generator power factor
 * type: real float
 * indexed
 */
#define GENERATOR_WPF "GENERATOR_WPF"

/**
 * Generator vset
 * type: real float
 * indexed
 */
#define GENERATOR_VSET "GENERATOR_VSET"

/**
 * Generator mq
 * type: real float
 * indexed
 */
#define GENERATOR_MQ "GENERATOR_MQ"

/**
 * Generator kpv
 * type: real float
 * indexed
 */
#define GENERATOR_KPV "GENERATOR_KPV"

/**
 * Generator kiv
 * type: real float
 * indexed
 */
#define GENERATOR_KIV "GENERATOR_KIV"

/**
 * Generator mp
 * type: real float
 * indexed
 */
#define GENERATOR_MP "GENERATOR_MP"

/**
 * Generator kppmax
 * type: real float
 * indexed
 */
#define GENERATOR_KPPMAX "GENERATOR_KPPMAX"

/**
 * Generator kipmax
 * type: real float
 * indexed
 */
#define GENERATOR_KIPMAX "GENERATOR_KIPMAX"

/**
 * Generator pmax
 * type: real float
 * indexed
 */
#define GENERATOR_PMAX "GENERATOR_PMAX"

/**
 * Generator pmin
 * type: real float
 * indexed
 */
#define GENERATOR_PMIN "GENERATOR_PMIN"

/**
 * Generator emax
 * type: real float
 * indexed
 */
#define GENERATOR_EMAX "GENERATOR_EMAX"

/**
 * Generator emin
 * type: real float
 * indexed
 */
#define GENERATOR_EMIN "GENERATOR_EMIN"

/**
 * Generator tpf
 * type: real float
 * indexed
 */
#define GENERATOR_TPF "GENERATOR_TPF"

/**
 * Generator imax
 * type: real float
 * indexed
 */
#define GENERATOR_IMAX "GENERATOR_IMAX"

/**
 * Generator qmax
 * type: real float
 * indexed
 */
#define GENERATOR_QMAX "GENERATOR_QMAX"

/**
 * Generator qmin
 * type: real float
 * indexed
 */
#define GENERATOR_QMIN "GENERATOR_QMIN"

/**
 * Generator kpqmax
 * type: real float
 * indexed
 */
#define GENERATOR_KPQMAX "GENERATOR_KPQMAX"

/**
 * Generator kiqmax
 * type: real float
 * indexed
 */
#define GENERATOR_KIQMAX "GENERATOR_KIQMAX"

/**
 * Generator tqf
 * type: real float
 * indexed
 */
#define GENERATOR_TQF "GENERATOR_TQF"

/**
 * Generator tvf
 * type: real float
 * indexed
 */
#define GENERATOR_TVF "GENERATOR_TVF"

/**
 * Generator vflag
 * type: int
 * indexed
 */
#define GENERATOR_VFLAG "GENERATOR_VFLAG"

// start generator REGCA parameters here
/**
 * Generator REGCA Lvplsw
 * type: integer
 * indexed
 */
#define GENERATOR_REGC_LVPLSW "GENERATOR_REGC_LVPLSW"

/**
 * Generator REGCA tg
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_TG "GENERATOR_REGC_TG"


/**
 * Generator REGCA Rrpwr
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_RRPWR "GENERATOR_REGC_RRPWR"

/**
 * Generator REGCA Brkpt
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_BRKPT "GENERATOR_REGC_BRKPT"

/**
 * Generator REGCA Zerox
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_ZEROX "GENERATOR_REGC_ZEROX"

/**
 * Generator REGCA Lvpl1
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_LVPL1 "GENERATOR_REGC_LVPL1"

/**
 * Generator REGCA Volim
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_VOLIM "GENERATOR_REGC_VOLIM"

/**
 * Generator REGCA Lvpnt1
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_LVPNT1 "GENERATOR_REGC_LVPNT1"

/**
 * Generator REGCA Lvpnt0
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_LVPNT0 "GENERATOR_REGC_LVPNT0"

/**
 * Generator REGCA lolim
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_LOLIM "GENERATOR_REGC_LOLIM"

/**
 * Generator REGCA Tfltr
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_TFLTR "GENERATOR_REGC_TFLTR"

/**
 * Generator REGCA Khv
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_KHV "GENERATOR_REGC_KHV"

/**
 * Generator REGCA iqrmax
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_IQRMAX "GENERATOR_REGC_IQRMAX"

/**
 * Generator REGCA iqrmin
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_IQRMIN "GENERATOR_REGC_IQRMIN"

/**
 * Generator REGCA accel
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_ACCEL "GENERATOR_REGC_ACCEL"

/**
 * Generator REGCB te
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_TE "GENERATOR_REGC_TE"

/**
 * Generator REGCB imax
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_IMAX "GENERATOR_REGC_IMAX"

/**
 * Generator REGCB RateFlag
 * type: integer
 * indexed
 */
#define GENERATOR_REGC_RATEFLAG "GENERATOR_REGC_RATEFLAG"

/**
 * Generator REGCB PQFLAG
 * type: integer
 * indexed
 */
#define GENERATOR_REGC_PQFLAG "GENERATOR_REGC_PQFLAG"

/**
 * Generator REGCC kip
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_KIP "GENERATOR_REGC_KIP"

/**
 * Generator REGCC kii
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_KII "GENERATOR_REGC_KII"

/**
 * Generator REGCC kipll
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_KIPLL "GENERATOR_REGC_KIPLL"

/**
 * Generator REGCC kppll
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_KPPLL "GENERATOR_REGC_KPPLL"

/**
 * Generator REGCC wmax
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_WMAX "GENERATOR_REGC_WMAX"

/**
 * Generator REGCC wmin
 * type: real float
 * indexed
 */
#define GENERATOR_REGC_WMIN "GENERATOR_REGC_WMIN"

#define HAS_PLANT_CONTROLLER "HAS_PLANT_CONTROLLER"

#define PLANT_CONTROLLER_MODEL "PLANT_CONTROLLER_MODEL"

// start generator REPCA parameters here

/**
 * Generator REPCA remote bus
 * type: integer
 * indexed
 */
#define GENERATOR_REPCA_IREG "GENERATOR_REPCA_IREG"

/**
 * Generator REPCA montiored branch FROM bus number for line compensation
 * type: integer
 * indexed
 */
#define GENERATOR_REPCA_BRCH_BUS_FROM "GENERATOR_REPCA_BRCH_BUS_FROM"

/**
 * Generator REPCA montiored branch TO bus number for line compensation
 * type: integer
 * indexed
 */
#define GENERATOR_REPCA_BRCH_BUS_TO "GENERATOR_REPCA_BRCH_BUS_TO"

/**
 * Generator REPCA montiored branch CKT for line compensation
 * type: string
 * indexed
 */
#define GENERATOR_REPCA_BRCH_CKT "GENERATOR_REPCA_BRCH_CKT"

/**
 * Generator REPCA VC flag (droop falg)
 * type: integer
 * indexed
 */
#define GENERATOR_REPCA_VC_FLAG "GENERATOR_REPCA_VC_FLAG"

/**
 * Generator REPCA Ref flag (flag for V or Q control)
 * type: integer
 * indexed
 */
#define GENERATOR_REPCA_REF_FLAG "GENERATOR_REPCA_REF_FLAG"

/**
 * Generator REPCA F flag to disable frequency control, 1: enable, 0:disable
 * type: integer
 * indexed
 */
#define GENERATOR_REPCA_F_FLAG "GENERATOR_REPCA_F_FLAG"

// finished REPCA M parameters definition, starts J parameters definition

/**
 * Generator REPCA Tfltr
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_TFLTR "GENERATOR_REPCA_TFLTR"

/**
 * Generator REPCA Kp
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_KP "GENERATOR_REPCA_KP"

/**
 * Generator REPCA Ki
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_KI "GENERATOR_REPCA_KI"

/**
 * Generator REPCA Tft
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_TFT "GENERATOR_REPCA_TFT"

/**
 * Generator REPCA Tfv
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_TFV "GENERATOR_REPCA_TFV"

/**
 * Generator REPCA Vfrz
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_VFRZ "GENERATOR_REPCA_VFRZ"

/**
 * Generator REPCA Rc
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_RC "GENERATOR_REPCA_RC"

/**
 * Generator REPCA Xc
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_XC "GENERATOR_REPCA_XC"

/**
 * Generator REPCA Kc
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_KC "GENERATOR_REPCA_KC"

/**
 * Generator REPCA emax
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_EMAX "GENERATOR_REPCA_EMAX"

/**
 * Generator REPCA emin
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_EMIN "GENERATOR_REPCA_EMIN"

/**
 * Generator REPCA dbd1
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_DBD1 "GENERATOR_REPCA_DBD1"

/**
 * Generator REPCA dbd2
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_DBD2 "GENERATOR_REPCA_DBD2"

/**
 * Generator REPCA qmax
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_QMAX "GENERATOR_REPCA_QMAX"

/**
 * Generator REPCA qmin
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_QMIN "GENERATOR_REPCA_QMIN"

/**
 * Generator REPCA kpg
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_KPG "GENERATOR_REPCA_KPG"

/**
 * Generator REPCA kig
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_KIG "GENERATOR_REPCA_KIG"

/**
 * Generator REPCA tp
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_TP "GENERATOR_REPCA_TP"

/**
 * Generator REPCA fdbd1
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_FDBD1 "GENERATOR_REPCA_FDBD1"

/**
 * Generator REPCA fdbd2
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_FDBD2 "GENERATOR_REPCA_FDBD2"

/**
 * Generator REPCA femax
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_FEMAX "GENERATOR_REPCA_FEMAX"

/**
 * Generator REPCA femin
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_FEMIN "GENERATOR_REPCA_FEMIN"

/**
 * Generator REPCA pmax
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_PMAX "GENERATOR_REPCA_PMAX"

/**
 * Generator REPCA pmin
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_PMIN "GENERATOR_REPCA_PMIN"

/**
 * Generator REPCA tg
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_TG "GENERATOR_REPCA_TG"

/**
 * Generator REPCA ddn
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_DDN "GENERATOR_REPCA_DDN"

/**
 * Generator REPCA dup
 * type: real float
 * indexed
 */
#define GENERATOR_REPCA_DUP "GENERATOR_REPCA_DUP"

// finished generator REPCA parameters definition


// start generator REECA parameters here

/**
 * Generator REECA REMOTE CONTROL BUS
 * type: integer
 * indexed
 */
#define GENERATOR_REECA_IREG "GENERATOR_REECA_IREG"

/**
 * Generator REECA PFFLAG
 * type: integer
 * indexed
 */
#define GENERATOR_REECA_PFFLAG "GENERATOR_REECA_PFFLAG"

/**
 * Generator REECA VFLAG
 * type: integer
 * indexed
 */
#define GENERATOR_REECA_VFLAG "GENERATOR_REECA_VFLAG"

/**
 * Generator REECA QFLAG
 * type: integer
 * indexed
 */
#define GENERATOR_REECA_QFLAG "GENERATOR_REECA_QFLAG"

/**
 * Generator REECA PFLAG
 * type: integer
 * indexed
 */
#define GENERATOR_REECA_PFLAG "GENERATOR_REECA_PFLAG"

/**
 * Generator REECA PQFLAG
 * type: integer
 * indexed
 */
#define GENERATOR_REECA_PQFLAG "GENERATOR_REECA_PQFLAG"

// finished M parameters definition for REECA, start J parameters defintion

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VDIP "GENERATOR_REECA_VDIP"


/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VUP "GENERATOR_REECA_VUP"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_TRV "GENERATOR_REECA_TRV"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_DBD1 "GENERATOR_REECA_DBD1"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_DBD2 "GENERATOR_REECA_DBD2"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_KQV "GENERATOR_REECA_KQV"


/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IQH1 "GENERATOR_REECA_IQH1"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IQL1 "GENERATOR_REECA_IQL1"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VREF0 "GENERATOR_REECA_VREF0"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IQFRZ "GENERATOR_REECA_IQFRZ"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_THLD "GENERATOR_REECA_THLD"


/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_THLD2 "GENERATOR_REECA_THLD2"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_TP "GENERATOR_REECA_TP"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_QMAX "GENERATOR_REECA_QMAX"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_QMIN "GENERATOR_REECA_QMIN"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VMAX "GENERATOR_REECA_VMAX"


/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VMIN "GENERATOR_REECA_VMIN"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_KQP "GENERATOR_REECA_KQP"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_KQI "GENERATOR_REECA_KQI"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_KVP "GENERATOR_REECA_KVP"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_KVI "GENERATOR_REECA_KVI"


/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VBIAS "GENERATOR_REECA_VBIAS"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_TIQ "GENERATOR_REECA_TIQ"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_DPMAX "GENERATOR_REECA_DPMAX"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_DPMIN "GENERATOR_REECA_DPMIN"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_PMAX "GENERATOR_REECA_PMAX"


/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_PMIN "GENERATOR_REECA_PMIN"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IMAX "GENERATOR_REECA_IMAX"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_TPORD "GENERATOR_REECA_TPORD"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VQ1 "GENERATOR_REECA_VQ1"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IQ1 "GENERATOR_REECA_IQ1"


/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VQ2 "GENERATOR_REECA_VQ2"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IQ2 "GENERATOR_REECA_IQ2"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VQ3 "GENERATOR_REECA_VQ3"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IQ3 "GENERATOR_REECA_IQ3"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VQ4 "GENERATOR_REECA_VQ4"


/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IQ4 "GENERATOR_REECA_IQ4"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VP1 "GENERATOR_REECA_VP1"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IP1 "GENERATOR_REECA_IP1"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VP2 "GENERATOR_REECA_VP2"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IP2 "GENERATOR_REECA_IP2"


/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VP3 "GENERATOR_REECA_VP3"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IP3 "GENERATOR_REECA_IP3"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_VP4 "GENERATOR_REECA_VP4"

/**
 * Generator REECA 
 * type: real float
 * indexed
 */
#define GENERATOR_REECA_IP4 "GENERATOR_REECA_IP4"

// finished generator REECA parameters definition

// Generator sequence parameters

/**
 * Generator positive sequence resistance
 * type: real float
 * indexed
 */
#define GENERATOR_SEQ_ZRPOS "GENERATOR_SEQ_ZRPOS"

/**
 * Generator positive sequence reactance
 * type: real float
 * indexed
 */
#define GENERATOR_SEQ_ZXPOS "GENERATOR_SEQ_ZXPOS"

/**
 * Generator negative sequence resistance
 * type: real float
 * indexed
 */
#define GENERATOR_SEQ_ZRNEG "GENERATOR_SEQ_ZRNEG"

/**
 * Generator negative sequence reactance
 * type: real float
 * indexed
 */
#define GENERATOR_SEQ_ZXNEG "GENERATOR_SEQ_ZXNEG"

/**
 * Generator zero sequence resistance
 * type: real float
 * indexed
 */
#define GENERATOR_SEQ_RZERO "GENERATOR_SEQ_RZERO"

/**
 * Generator zero sequence reactance
 * type: real float
 * indexed
 */
#define GENERATOR_SEQ_XZERO "GENERATOR_SEQ_XZERO"

// GENERATOR COST PARAMETERS
/**
 * Generator cost model type:
 *  1 piecewise linear
 *  2 polynomial
 *  type: integer
 *  indexed
 */
#define GENCOST_MODEL_TYPE "GENCOST_MODEL_TYPE"

/**
 * Startup cost in US dollars
 * type: real float
 * indexed
 */
#define GENCOST_STARTUP "GENCOST_STARTUP"

/**
 * Shutdown cost in US dollars
 * type: real float
 * indexed
 */
#define GENCOST_SHUTDOWN "GENCOST_SHUTDOWN"

/**
 * Number of parameters in model
 * type: integer
 * indexed
 */
#define GENCOST_NUM_PARAMS "GENCOST_NUM_PARAMS"

/**
 * Because the maximum number of parameters is unknown, it is not possible to
 * predefine these variables. However, they have the format
 * Model = 1
 * p0, p1, ..., pN
 * f0, f1, ..., fN
 * N is the number parameters in the model, pi signifies the
 * break/endpoints of the intervals and fi are the corresponding values at pi.
 * Model = 2
 * cN, ..., c1, c0
 * The ci are the coefficients of a polynomial of order N
 * cN*p**N + ... + c1*p + c0
 * The parameters are denoted by the labels
 * GENCOST_PARAM_PN
 * GENCOST_PARAM_FN
 * GENCOST_PARAN_CN
 *
 * The corresponding cost parameters for reactive power, if included are
 * GENCOST_PARAM_R_PN
 * GENCOST_PARAM_R_FN
 * GENCOST_PARAN_R_CN
 */

/**
 * Non-blank alphanumeric string to identify model being used for pss	
 * type: string
 * indexed
 */
#define PSSSIM_MODEL "PSSSIM_MODEL"

/**
 * Non-blank alphanumeric string to identify model being used for generator	
 * type: string
 * indexed
 */
#define GENERATOR_MODEL "GENERATOR_MODEL"

/**
 * Nominal power of the generator
 * type: float
 * indexed
 */
#define GENERATOR_NOM_POWER_CURRENT "GENERATOR_NOM_POWER_CURRENT"

#endif /* _GENERATOR_VAR_HPP_ */
