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
 *  multiple generators can be associated with a bus and multiple transmission
 *  elements can be associated with a branch. The variables that have an associated
 *  index are denoted with the keyword "indexed".
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

// LOAD DATA

/**
 * Number of loads on bus
 * type: integer
 */
#define LOAD_NUMBER "LOAD_NUMBER"

/**
 * The Bus number to which the load is connected
 * type: integer
 * indexed
 */
#define LOAD_BUSNUMBER "LOAD_BUSNUMBER"

/**
 * Non-blank alphanumeric identifier to distinguish different loads connected to the same bus.
 * Default value: ’1’
 * type: string
 * indexed
 */
#define LOAD_ID "LOAD_ID"

/**
 * Load status
 *  1: in-service
 *  0: out-of-service
 * Default value is 1
 * type: integer
 * indexed
 */
#define LOAD_STATUS "LOAD_STATUS"

/**
 * Area to which the load is assigned
 * type: integer
 * indexed
 */
#define LOAD_AREA "LOAD_AREA"

/**
 * Zone to which the load is assigned
 * type: integer
 * indexed
 */
#define LOAD_ZONE "LOAD_ZONE"

/**
 * Active power component of constant power load; entered in MW
 * type: real float
 * indexed
 */
#define LOAD_PL "LOAD_PL"

/**
 * Reactive power component of constant power load; entered in MVar
 * type: real float
 * indexed
 */
#define LOAD_QL "LOAD_QL"

/**
 * Active power component of constant current load; entered in MW 
 * type: real float
 * indexed
 */
#define LOAD_IP "LOAD_IP"

/**
 * Reactive power component of constant current load; entered in Mvar
 * type: real float
 * indexed
 */
#define LOAD_IQ "LOAD_IQ"

/**
 * Active power component of constant admittance load; entered in MW
 * type: real float
 * indexed
 */
#define LOAD_YP "LOAD_YP"

/**
 * Reactive power component of constant admittance load; entered in MVar
 * type: real float
 * indexed
 */
#define LOAD_YQ "LOAD_YQ"

/**
 * Owner to which load is assigned
 * type: integer
 * indexed
 */
#define LOAD_OWNER "LOAD_OWNER"

/**
 * Load scaling flag
 * type: integer
 * indexed
 */
#define LOAD_SCALE "LOAD_SCALE"

/**
 * Interruptible load flag
 * type: integer
 * indexed
 */
#define LOAD_INTRPT "LOAD_INTRPT"

/**
 * Alphanumeric string describing load model
 * type: string
 * indexed
 */
#define LOAD_MODEL "LOAD_MODEL"

/**
 * Load IT
 * type: integer
 * indexed
 */
#define LOAD_IT "LOAD_IT"

/**
 * Load RA
 * type: float
 * indexed
 */
#define LOAD_RA "LOAD_RA"

/**
 * Load XA
 * type: float
 * indexed
 */
#define LOAD_XA "LOAD_XA"

/**
 * Load XM
 * type: float
 * indexed
 */
#define LOAD_XM "LOAD_XM"

/**
 * Load R1
 * type: float
 * indexed
 */
#define LOAD_R1 "LOAD_R1"

/**
 * Load X1
 * type: float
 * indexed
 */
#define LOAD_X1 "LOAD_X1"

/**
 * Load R2
 * type: float
 * indexed
 */
#define LOAD_R2 "LOAD_R2"

/**
 * Load X1
 * type: float
 * indexed
 */
#define LOAD_X2 "LOAD_X2"

/**
 * Load E1
 * type: float
 * indexed
 */
#define LOAD_E1 "LOAD_E1"

/**
 * Load SE1
 * type: float
 * indexed
 */
#define LOAD_SE1 "LOAD_SE1"

/**
 * Load E2
 * type: float
 * indexed
 */
#define LOAD_E2 "LOAD_E2"

/**
 * Load SE2
 * type: float
 * indexed
 */
#define LOAD_SE2 "LOAD_SE2"

/**
 * Load MBASE
 * type: float
 * indexed
 */
#define LOAD_MBASE "LOAD_MBASE"

/**
 * Load PMULT
 * type: float
 * indexed
 */
#define LOAD_PMULT "LOAD_PMULT"

/**
 * Load H
 * type: float
 * indexed
 */
#define LOAD_H "LOAD_H"

/**
 * Load VI
 * type: float
 * indexed
 */
#define LOAD_VI "LOAD_VI"

/**
 * Load TI
 * type: float
 * indexed
 */
#define LOAD_TI "LOAD_TI"

/**
 * Load TB
 * type: float
 * indexed
 */
#define LOAD_TB "LOAD_TB"

/**
 * Load A
 * type: float
 * indexed
 */
#define LOAD_A "LOAD_A"

/**
 * Load B
 * type: float
 * indexed
 */
#define LOAD_B "LOAD_B"

/**
 * Load D
 * type: float
 * indexed
 */
#define LOAD_D "LOAD_D"

/**
 * Load E
 * type: float
 * indexed
 */
#define LOAD_E "LOAD_E"

/**
 * Load C0
 * type: float
 * indexed
 */
#define LOAD_C0 "LOAD_C0"

/**
 * Load TNOM
 * type: float
 * indexed
 */
#define LOAD_TNOM "LOAD_TNOM"

/**
 * Load TSTALL
 * type: float
 * indexed
 */
#define LOAD_TSTALL "LOAD_TSTALL"

/**
 * Load TRESTART
 * type: float
 * indexed
 */
#define LOAD_TRESTART "LOAD_TRESTART"

/**
 * Load TV
 * type: float
 * indexed
 */
#define LOAD_TV "LOAD_TV"

/**
 * Load TF
 * type: float
 * indexed
 */
#define LOAD_TF "LOAD_TF"

/**
 * Load COMPLF
 * type: float
 * indexed
 */
#define LOAD_COMPLF "LOAD_COMPLF"

/**
 * Load COMPPF
 * type: float
 * indexed
 */
#define LOAD_COMPPF "LOAD_COMPPF"

/**
 * Load VSTALL
 * type: float
 * indexed
 */
#define LOAD_VSTALL "LOAD_VSTALL"

/**
 * Load RSTALL
 * type: float
 * indexed
 */
#define LOAD_RSTALL "LOAD_RSTALL"

/**
 * Load XSTALL
 * type: float
 * indexed
 */
#define LOAD_XSTALL "LOAD_XSTALL"

/**
 * Load LFADJ
 * type: float
 * indexed
 */
#define LOAD_LFADJ "LOAD_LFADJ"

/**
 * Load KP1
 * type: float
 * indexed
 */
#define LOAD_KP1 "LOAD_KP1"

/**
 * Load NP1
 * type: float
 * indexed
 */
#define LOAD_NP1 "LOAD_NP1"

/**
 * Load KQ1
 * type: float
 * indexed
 */
#define LOAD_KQ1 "LOAD_KQ1"

/**
 * Load NQ1
 * type: float
 * indexed
 */
#define LOAD_NQ1 "LOAD_NQ1"

/**
 * Load KP2
 * type: float
 * indexed
 */
#define LOAD_KP2 "LOAD_KP2"

/**
 * Load NP2
 * type: float
 * indexed
 */
#define LOAD_NP2 "LOAD_NP2"

/**
 * Load KQ2
 * type: float
 * indexed
 */
#define LOAD_KQ2 "LOAD_KQ2"

/**
 * Load NQ2
 * type: float
 * indexed
 */
#define LOAD_NQ2 "LOAD_NQ2"

/**
 * Load VBRK
 * type: float
 * indexed
 */
#define LOAD_VBRK "LOAD_VBRK"

/**
 * Load FRST
 * type: float
 * indexed
 */
#define LOAD_FRST "LOAD_FRST"

/**
 * Load VRST
 * type: float
 * indexed
 */
#define LOAD_VRST "LOAD_VRST"

/**
 * Load CMPKPF
 * type: float
 * indexed
 */
#define LOAD_CMPKPF "LOAD_CMPKPF"

/**
 * Load CMPKQF
 * type: float
 * indexed
 */
#define LOAD_CMPKQF "LOAD_CMPKQF"

/**
 * Load VC1OFF
 * type: float
 * indexed
 */
#define LOAD_VC1OFF "LOAD_VC1OFF"

/**
 * Load VC2OFF
 * type: float
 * indexed
 */
#define LOAD_VC2OFF "LOAD_VC2OFF"

/**
 * Load VC1ON
 * type: float
 * indexed
 */
#define LOAD_VC1ON "LOAD_VC1ON"

/**
 * Load VC2ON
 * type: float
 * indexed
 */
#define LOAD_VC2ON "LOAD_VC2ON"

/**
 * Load TTH
 * type: float
 * indexed
 */
#define LOAD_TTH "LOAD_TTH"

/**
 * Load TH1T
 * type: float
 * indexed
 */
#define LOAD_TH1T "LOAD_TH1T"

/**
 * Load TH2T
 * type: float
 * indexed
 */
#define LOAD_TH2T "LOAD_TH2T"

/**
 * Load FUVR
 * type: float
 * indexed
 */
#define LOAD_FUVR "LOAD_FUVR"

/**
 * Load UVTR1
 * type: float
 * indexed
 */
#define LOAD_UVTR1 "LOAD_UVTR1"

/**
 * Load TTR1
 * type: float
 * indexed
 */
#define LOAD_TTR1 "LOAD_TTR1"

/**
 * Load UVTR2
 * type: float
 * indexed
 */
#define LOAD_UVTR2 "LOAD_UVTR2"

/**
 * Load TTR2
 * type: float
 * indexed
 */
#define LOAD_TTR2 "LOAD_TTR2"

/**
 * Load A1
 * type: float
 * indexed
 */
#define LOAD_A1 "LOAD_A1"

/**
 * Load A2
 * type: float
 * indexed
 */
#define LOAD_A2 "LOAD_A2"

/**
 * Load A3
 * type: float
 * indexed
 */
#define LOAD_A3 "LOAD_A3"

/**
 * Load A4
 * type: float
 * indexed
 */
#define LOAD_A4 "LOAD_A4"

/**
 * Load A5
 * type: float
 * indexed
 */
#define LOAD_A5 "LOAD_A5"

/**
 * Load A6
 * type: float
 * indexed
 */
#define LOAD_A6 "LOAD_A6"

/**
 * Load A7
 * type: float
 * indexed
 */
#define LOAD_A7 "LOAD_A7"

/**
 * Load A8
 * type: float
 * indexed
 */
#define LOAD_A8 "LOAD_A8"

/**
 * Load N1
 * type: float
 * indexed
 */
#define LOAD_N1 "LOAD_N1"

/**
 * Load N2
 * type: float
 * indexed
 */
#define LOAD_N2 "LOAD_N2"

/**
 * Load N3
 * type: float
 * indexed
 */
#define LOAD_N3 "LOAD_N3"

/**
 * Load N4
 * type: float
 * indexed
 */
#define LOAD_N4 "LOAD_N4"

/**
 * Load N5
 * type: float
 * indexed
 */
#define LOAD_N5 "LOAD_N5"

/**
 * Load N6
 * type: float
 * indexed
 */
#define LOAD_N6 "LOAD_N6"

/**
 * Load MVA
 * type: float
 * indexed
 */
#define LOAD_MVA "LOAD_MVA"

/**
 * Load BSS
 * type: float
 * indexed
 */
#define LOAD_BSS "LOAD_BSS"

/**
 * Load RFDR
 * type: float
 * indexed
 */
#define LOAD_RFDR "LOAD_RFDR"

/**
 * Load XFDR
 * type: float
 * indexed
 */
#define LOAD_XFDR "LOAD_XFDR"

/**
 * Load FB
 * type: float
 * indexed
 */
#define LOAD_FB "LOAD_FB"

/**
 * Load XXF
 * type: float
 * indexed
 */
#define LOAD_XXF "LOAD_XXF"

/**
 * Load TFIXHS
 * type: float
 * indexed
 */
#define LOAD_TFIXHS "LOAD_TFIXHS"

/**
 * Load TFIXLS
 * type: float
 * indexed
 */
#define LOAD_TFIXLS "LOAD_TFIXLS"

/**
 * Load LTC
 * type: float
 * indexed
 */
#define LOAD_LTC "LOAD_LTC"

/**
 * Load TMIN
 * type: float
 * indexed
 */
#define LOAD_TMIN "LOAD_TMIN"

/**
 * Load TMAX
 * type: float
 * indexed
 */
#define LOAD_TMAX "LOAD_TMAX"

/**
 * Load STEP
 * type: float
 * indexed
 */
#define LOAD_STEP "LOAD_STEP"

/**
 * Load VMIN
 * type: float
 * indexed
 */
#define LOAD_VMIN "LOAD_VMIN"

/**
 * Load VMAX
 * type: float
 * indexed
 */
#define LOAD_VMAX "LOAD_VMAX"

/**
 * Load TDEL
 * type: float
 * indexed
 */
#define LOAD_TDEL "LOAD_TDEL"

/**
 * Load TTAP
 * type: float
 * indexed
 */
#define LOAD_TTAP "LOAD_TTAP"

/**
 * Load RCOMP
 * type: float
 * indexed
 */
#define LOAD_RCOMP "LOAD_RCOMP"

/**
 * Load XCOMP
 * type: float
 * indexed
 */
#define LOAD_XCOMP "LOAD_XCOMP"

/**
 * Load FMA
 * type: float
 * indexed
 */
#define LOAD_FMA "LOAD_FMA"

/**
 * Load FMB
 * type: float
 * indexed
 */
#define LOAD_FMB "LOAD_FMB"

/**
 * Load FMC
 * type: float
 * indexed
 */
#define LOAD_FMC "LOAD_FMC"

/**
 * Load FMD
 * type: float
 * indexed
 */
#define LOAD_FMD "LOAD_FMD"

/**
 * Load FEL
 * type: float
 * indexed
 */
#define LOAD_FEL "LOAD_FEL"

/**
 * Load PFEL
 * type: float
 * indexed
 */
#define LOAD_PFEL "LOAD_PFEL"

/**
 * Load VD1
 * type: float
 * indexed
 */
#define LOAD_VD1 "LOAD_VD1"

/**
 * Load VD2
 * type: float
 * indexed
 */
#define LOAD_VD2 "LOAD_VD2"

/**
 * Load FRCEL
 * type: float
 * indexed
 */
#define LOAD_FRCEL "LOAD_FRCEL"

/**
 * Load PFS
 * type: float
 * indexed
 */
#define LOAD_PFS "LOAD_PFS"

/**
 * Load P1E
 * type: float
 * indexed
 */
#define LOAD_P1E "LOAD_P1E"

/**
 * Load P1C
 * type: float
 * indexed
 */
#define LOAD_P1C "LOAD_P1C"

/**
 * Load P2E
 * type: float
 * indexed
 */
#define LOAD_P2E "LOAD_P2E"

/**
 * Load P2C
 * type: float
 * indexed
 */
#define LOAD_P2C "LOAD_P2C"

/**
 * Load PFREQ
 * type: float
 * indexed
 */
#define LOAD_PFREQ "LOAD_PFREQ"

/**
 * Load Q1E
 * type: float
 * indexed
 */
#define LOAD_Q1E "LOAD_Q1E"

/**
 * Load Q1C
 * type: float
 * indexed
 */
#define LOAD_Q1C "LOAD_Q1C"

/**
 * Load Q2E
 * type: float
 * indexed
 */
#define LOAD_Q2E "LOAD_Q2E"

/**
 * Load Q2C
 * type: float
 * indexed
 */
#define LOAD_Q2C "LOAD_Q2C"

/**
 * Load QFREQ
 * type: float
 * indexed
 */
#define LOAD_QFREQ "LOAD_QFREQ"

/**
 * Load MTPA
 * type: integer
 * indexed
 */
#define LOAD_MTPA "LOAD_MTPA"

/**
 * Load LFMA
 * type: float
 * indexed
 */
#define LOAD_LFMA "LOAD_LFMA"

/**
 * Load RSA
 * type: float
 * indexed
 */
#define LOAD_RSA "LOAD_RSA"

/**
 * Load LSA
 * type: float
 * indexed
 */
#define LOAD_LSA "LOAD_LSA"

/**
 * Load LPA
 * type: float
 * indexed
 */
#define LOAD_LPA "LOAD_LPA"

/**
 * Load LPPA
 * type: float
 * indexed
 */
#define LOAD_LPPA "LOAD_LPPA"

/**
 * Load TPOA
 * type: float
 * indexed
 */
#define LOAD_TPOA "LOAD_TPOA"

/**
 * Load TPPOA
 * type: float
 * indexed
 */
#define LOAD_TPPOA "LOAD_TPPOA"

/**
 * Load HA
 * type: float
 * indexed
 */
#define LOAD_HA "LOAD_HA"

/**
 * Load ETRQA
 * type: float
 * indexed
 */
#define LOAD_ETRQA "LOAD_ETRQA"

/**
 * Load VTR1A
 * type: float
 * indexed
 */
#define LOAD_VTR1A "LOAD_VTR1A"

/**
 * Load TTR1A
 * type: float
 * indexed
 */
#define LOAD_TTR1A "LOAD_TTR1A"

/**
 * Load FTR1A
 * type: float
 * indexed
 */
#define LOAD_FTR1A "LOAD_FTR1A"

/**
 * Load VRC1A
 * type: float
 * indexed
 */
#define LOAD_VRC1A "LOAD_VRC1A"

/**
 * Load TRC1A
 * type: float
 * indexed
 */
#define LOAD_TRC1A "LOAD_TRC1A"

/**
 * Load VTR2A
 * type: float
 * indexed
 */
#define LOAD_VTR2A "LOAD_VTR2A"

/**
 * Load TTR2A
 * type: float
 * indexed
 */
#define LOAD_TTR2A "LOAD_TTR2A"

/**
 * Load FTR2A
 * type: float
 * indexed
 */
#define LOAD_FTR2A "LOAD_FTR2A"

/**
 * Load VRC2A
 * type: float
 * indexed
 */
#define LOAD_VRC2A "LOAD_VRC2A"

/**
 * Load TRC2A
 * type: float
 * indexed
 */
#define LOAD_TRC2A "LOAD_TRC2A"

/**
 * Load MTPB
 * type: integer
 * indexed
 */
#define LOAD_MTPB "LOAD_MTPB"

/**
 * Load LFMB
 * type: float
 * indexed
 */
#define LOAD_LFMB "LOAD_LFMB"

/**
 * Load RSB
 * type: float
 * indexed
 */
#define LOAD_RSB "LOAD_RSB"

/**
 * Load LSB
 * type: float
 * indexed
 */
#define LOAD_LSB "LOAD_LSB"

/**
 * Load LPB
 * type: float
 * indexed
 */
#define LOAD_LPB "LOAD_LPB"

/**
 * Load LPPB
 * type: float
 * indexed
 */
#define LOAD_LPPB "LOAD_LPPB"

/**
 * Load TPOB
 * type: float
 * indexed
 */
#define LOAD_TPOB "LOAD_TPOB"

/**
 * Load TPPOB
 * type: float
 * indexed
 */
#define LOAD_TPPOB "LOAD_TPPOB"

/**
 * Load HB
 * type: float
 * indexed
 */
#define LOAD_HB "LOAD_HB"

/**
 * Load ETRQB
 * type: float
 * indexed
 */
#define LOAD_ETRQB "LOAD_ETRQB"

/**
 * Load VTR1B
 * type: float
 * indexed
 */
#define LOAD_VTR1B "LOAD_VTR1B"

/**
 * Load TTR1B
 * type: float
 * indexed
 */
#define LOAD_TTR1B "LOAD_TTR1B"

/**
 * Load FTR1B
 * type: float
 * indexed
 */
#define LOAD_FTR1B "LOAD_FTR1B"

/**
 * Load VRC1B
 * type: float
 * indexed
 */
#define LOAD_VRC1B "LOAD_VRC1B"

/**
 * Load TRC1B
 * type: float
 * indexed
 */
#define LOAD_TRC1B "LOAD_TRC1B"

/**
 * Load VTR2B
 * type: float
 * indexed
 */
#define LOAD_VTR2B "LOAD_VTR2B"

/**
 * Load TTR2B
 * type: float
 * indexed
 */
#define LOAD_TTR2B "LOAD_TTR2B"

/**
 * Load FTR2B
 * type: float
 * indexed
 */
#define LOAD_FTR2B "LOAD_FTR2B"

/**
 * Load VRC2B
 * type: float
 * indexed
 */
#define LOAD_VRC2B "LOAD_VRC2B"

/**
 * Load TRC2B
 * type: float
 * indexed
 */
#define LOAD_TRC2B "LOAD_TRC2B"

/**
 * Load MTPC
 * type: integer
 * indexed
 */
#define LOAD_MTPC "LOAD_MTPC"

/**
 * Load LFMC
 * type: float
 * indexed
 */
#define LOAD_LFMC "LOAD_LFMC"

/**
 * Load RSC
 * type: float
 * indexed
 */
#define LOAD_RSC "LOAD_RSC"

/**
 * Load LSC
 * type: float
 * indexed
 */
#define LOAD_LSC "LOAD_LSC"

/**
 * Load LPC
 * type: float
 * indexed
 */
#define LOAD_LPC "LOAD_LPC"

/**
 * Load LPPC
 * type: float
 * indexed
 */
#define LOAD_LPPC "LOAD_LPPC"

/**
 * Load TPOC
 * type: float
 * indexed
 */
#define LOAD_TPOC "LOAD_TPOC"

/**
 * Load TPPOC
 * type: float
 * indexed
 */
#define LOAD_TPPOC "LOAD_TPPOC"

/**
 * Load HC
 * type: float
 * indexed
 */
#define LOAD_HC "LOAD_HC"

/**
 * Load ETRQC
 * type: float
 * indexed
 */
#define LOAD_ETRQC "LOAD_ETRQC"

/**
 * Load VTR1C
 * type: float
 * indexed
 */
#define LOAD_VTR1C "LOAD_VTR1C"

/**
 * Load TTR1C
 * type: float
 * indexed
 */
#define LOAD_TTR1C "LOAD_TTR1C"

/**
 * Load FTR1C
 * type: float
 * indexed
 */
#define LOAD_FTR1C "LOAD_FTR1C"

/**
 * Load VRC1C
 * type: float
 * indexed
 */
#define LOAD_VRC1C "LOAD_VRC1C"

/**
 * Load TRC1C
 * type: float
 * indexed
 */
#define LOAD_TRC1C "LOAD_TRC1C"

/**
 * Load VTR2C
 * type: float
 * indexed
 */
#define LOAD_VTR2C "LOAD_VTR2C"

/**
 * Load TTR2C
 * type: float
 * indexed
 */
#define LOAD_TTR2C "LOAD_TTR2C"

/**
 * Load FTR2C
 * type: float
 * indexed
 */
#define LOAD_FTR2C "LOAD_FTR2C"

/**
 * Load VRC2C
 * type: float
 * indexed
 */
#define LOAD_VRC2C "LOAD_VRC2C"

/**
 * Load TRC2C
 * type: float
 * indexed
 */
#define LOAD_TRC2C "LOAD_TRC2C"

/**
 * Load MTPD
 * type: integer
 * indexed
 */
#define LOAD_MTPD "LOAD_MTPD"

/**
 * Load LFMD
 * type: float
 * indexed
 */
#define LOAD_LFMD "LOAD_LFMD"

/**
 * Load TRST
 * type: float
 * indexed
 */
#define LOAD_TRST "LOAD_TRST"

/**
 * Load VTR1
 * type: float
 * indexed
 */
#define LOAD_VTR1 "LOAD_VTR1"

/**
 * Load VTR2
 * type: float
 * indexed
 */
#define LOAD_VTR2 "LOAD_VTR2"

/**
 * Load LFM
 * type: float
 * indexed
 */
#define LOAD_LFM "LOAD_LFM"

/**
 * Load RS
 * type: float
 * indexed
 */
#define LOAD_RS "LOAD_RS"

/**
 * Load LS
 * type: float
 * indexed
 */
#define LOAD_LS "LOAD_LS"

/**
 * Load LP
 * type: float
 * indexed
 */
#define LOAD_LP "LOAD_LP"

/**
 * Load LPP
 * type: float
 * indexed
 */
#define LOAD_LPP "LOAD_LPP"

/**
 * Load TPO
 * type: float
 * indexed
 */
#define LOAD_TPO "LOAD_TPO"

/**
 * Load TPPO
 * type: float
 * indexed
 */
#define LOAD_TPPO "LOAD_TPPO"

/**
 * Load ETRQ
 * type: float
 * indexed
 */
#define LOAD_ETRQ "LOAD_ETRQ"

/**
 * Load FTR1
 * type: float
 * indexed
 */
#define LOAD_FTR1 "LOAD_FTR1"

/**
 * Load VRC1
 * type: float
 * indexed
 */
#define LOAD_VRC1 "LOAD_VRC1"

/**
 * Load TRC1
 * type: float
 * indexed
 */
#define LOAD_TRC1 "LOAD_TRC1"

/**
 * Load FTR2
 * type: float
 * indexed
 */
#define LOAD_FTR2 "LOAD_FTR2"

/**
 * Load VRC2
 * type: float
 * indexed
 */
#define LOAD_VRC2 "LOAD_VRC2"

/**
 * Load TRC2
 * type: float
 * indexed
 */
#define LOAD_TRC2 "LOAD_TRC2"

/**
 * Load MTP
 * type: integer
 * indexed
 */
#define LOAD_MTP "LOAD_MTP"

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

// GOVERNOR DATA
/**
 * Flag to indicate that governor is present
 * type: boolean
 * indexed
 */
#define HAS_GOVERNOR "HAS_GOVERNOR"

/**
 * Governor model
 * type: string
 * indexed
 */
#define GOVERNOR_MODEL "GOVERNOR_MODEL"

/**
 * Generator JBUS
 * type: integer
 * indexed
 */
#define GOVERNOR_JBUS "GOVERNOR_JBUS"

/**
 * Governor M
 * type: integer
 * indexed
 */
#define GOVERNOR_M "GOVERNOR_M"

/**
 * Governor K
 * type: real float
 * indexed
 */
#define GOVERNOR_K "GOVERNOR_K"

/**
 * Governor T1
 * type: real float
 * indexed
 */
#define GOVERNOR_T1 "GOVERNOR_T1"

/**
 * Governor T2
 * type: real float
 * indexed
 */
#define GOVERNOR_T2 "GOVERNOR_T2"

/**
 * Governor T3
 * type: real float
 * indexed
 */
#define GOVERNOR_T3 "GOVERNOR_T3"

/**
 * Governor UO
 * type: real float
 * indexed
 */
#define GOVERNOR_UO "GOVERNOR_UO"

/**
 * Governor UC
 * type: real float
 * indexed
 */
#define GOVERNOR_UC "GOVERNOR_UC"

/**
 * Governor PMAX
 * type: real float
 * indexed
 */
#define GOVERNOR_PMAX "GOVERNOR_PMAX"

/**
 * Governor PMIN
 * type: real float
 * indexed
 */
#define GOVERNOR_PMIN "GOVERNOR_PMIN"

/**
 * Governor T4
 * type: real float
 * indexed
 */
#define GOVERNOR_T4 "GOVERNOR_T4"

/**
 * Governor K1
 * type: real float
 * indexed
 */
#define GOVERNOR_K1 "GOVERNOR_K1"

/**
 * Governor K2
 * type: real float
 * indexed
 */
#define GOVERNOR_K2 "GOVERNOR_K2"

/**
 * Governor T5
 * type: real float
 * indexed
 */
#define GOVERNOR_T5 "GOVERNOR_T5"

/**
 * Governor K3
 * type: real float
 * indexed
 */
#define GOVERNOR_K3 "GOVERNOR_K3"

/**
 * Governor K4
 * type: real float
 * indexed
 */
#define GOVERNOR_K4 "GOVERNOR_K4"

/**
 * Governor T6
 * type: real float
 * indexed
 */
#define GOVERNOR_T6 "GOVERNOR_T6"

/**
 * Governor K5
 * type: real float
 * indexed
 */
#define GOVERNOR_K5 "GOVERNOR_K5"

/**
 * Governor K6
 * type: real float
 * indexed
 */
#define GOVERNOR_K6 "GOVERNOR_K6"

/**
 * Governor T7
 * type: real float
 * indexed
 */
#define GOVERNOR_T7 "GOVERNOR_T7"

/**
 * Governor K7
 * type: real float
 * indexed
 */
#define GOVERNOR_K7 "GOVERNOR_K7"

/**
 * Governor K8
 * type: real float
 * indexed
 */
#define GOVERNOR_K8 "GOVERNOR_K8"

/**
 * Governor DB1
 * type: real float
 * indexed
 */
#define GOVERNOR_DB1 "GOVERNOR_DB1"

/**
 * Governor ERR
 * type: real float
 * indexed
 */
#define GOVERNOR_ERR "GOVERNOR_ERR"

/**
 * Governor DB2
 * type: real float
 * indexed
 */
#define GOVERNOR_DB2 "GOVERNOR_DB2"

/**
 * Governor GV1
 * type: real float
 * indexed
 */
#define GOVERNOR_GV1 "GOVERNOR_GV1"

/**
 * Governor PGV1
 * type: real float
 * indexed
 */
#define GOVERNOR_PGV1 "GOVERNOR_PGV1"

/**
 * Governor GV2
 * type: real float
 * indexed
 */
#define GOVERNOR_GV2 "GOVERNOR_GV2"

/**
 * Governor PGV2
 * type: real float
 * indexed
 */
#define GOVERNOR_PGV2 "GOVERNOR_PGV2"

/**
 * Governor GV3
 * type: real float
 * indexed
 */
#define GOVERNOR_GV3 "GOVERNOR_GV3"

/**
 * Governor PGV3
 * type: real float
 * indexed
 */
#define GOVERNOR_PGV3 "GOVERNOR_PGV3"

/**
 * Governor GV4
 * type: real float
 * indexed
 */
#define GOVERNOR_GV4 "GOVERNOR_GV4"

/**
 * Governor PGV4
 * type: real float
 * indexed
 */
#define GOVERNOR_PGV4 "GOVERNOR_PGV4"

/**
 * Governor GV5
 * type: real float
 * indexed
 */
#define GOVERNOR_GV5 "GOVERNOR_GV5"

/**
 * Governor PGV5
 * type: real float
 * indexed
 */
#define GOVERNOR_PGV5 "GOVERNOR_PGV5"

/**
 * Governor IBLOCK
 * type: integer
 * indexed
 */
#define GOVERNOR_IBLOCK "GOVERNOR_IBLOCK"

/**
 * Governor RSELECT
 * type: real float
 * indexed
 */
#define GOVERNOR_RSELECT "GOVERNOR_RSELECT"

/**
 * Governor FLAGSWITCH
 * type: real float
 * indexed
 */
#define GOVERNOR_FLAGSWITCH "GOVERNOR_FLAGSWITCH"

/**
 * Governor R
 * type: real float
 * indexed
 */
#define GOVERNOR_R "GOVERNOR_R"

/**
 * Governor TPELEC
 * type: real float
 * indexed
 */
#define GOVERNOR_TPELEC "GOVERNOR_TPELEC"

/**
 * Governor MAXERR
 * type: real float
 * indexed
 */
#define GOVERNOR_MAXERR "GOVERNOR_MAXERR"

/**
 * Governor MINERR
 * type: real float
 * indexed
 */
#define GOVERNOR_MINERR "GOVERNOR_MINERR"

/**
 * Governor KPGOV
 * type: real float
 * indexed
 */
#define GOVERNOR_KPGOV "GOVERNOR_KPGOV"

/**
 * Governor KIGOV
 * type: real float
 * indexed
 */
#define GOVERNOR_KIGOV "GOVERNOR_KIGOV"

/**
 * Governor KDGOV
 * type: real float
 * indexed
 */
#define GOVERNOR_KDGOV "GOVERNOR_KDGOV"

/**
 * Governor TDGOV
 * type: real float
 * indexed
 */
#define GOVERNOR_TDGOV "GOVERNOR_TDGOV"

/**
 * Governor VMAX
 * type: real float
 * indexed
 */
#define GOVERNOR_VMAX "GOVERNOR_VMAX"

/**
 * Governor VMIN
 * type: real float
 * indexed
 */
#define GOVERNOR_VMIN "GOVERNOR_VMIN"

/**
 * Governor TACT
 * type: real float
 * indexed
 */
#define GOVERNOR_TACT "GOVERNOR_TACT"

/**
 * Governor KTURB
 * type: real float
 * indexed
 */
#define GOVERNOR_KTURB "GOVERNOR_KTURB"

/**
 * Governor WFNL
 * type: real float
 * indexed
 */
#define GOVERNOR_WFNL "GOVERNOR_WFNL"

/**
 * Governor TB
 * type: real float
 * indexed
 */
#define GOVERNOR_TB "GOVERNOR_TB"

/**
 * Governor TC
 * type: real float
 * indexed
 */
#define GOVERNOR_TC "GOVERNOR_TC"

/**
 * Governor TENG
 * type: real float
 * indexed
 */
#define GOVERNOR_TENG "GOVERNOR_TENG"

/**
 * Governor TFLOAD
 * type: real float
 * indexed
 */
#define GOVERNOR_TFLOAD "GOVERNOR_TFLOAD"

/**
 * Governor KPLOAD
 * type: real float
 * indexed
 */
#define GOVERNOR_KPLOAD "GOVERNOR_KPLOAD"

/**
 * Governor KILOAD
 * type: real float
 * indexed
 */
#define GOVERNOR_KILOAD "GOVERNOR_KILOAD"

/**
 * Governor LDREF
 * type: real float
 * indexed
 */
#define GOVERNOR_LDREF "GOVERNOR_LDREF"

/**
 * Governor DM
 * type: real float
 * indexed
 */
#define GOVERNOR_DM "GOVERNOR_DM"

/**
 * Governor ROPEN
 * type: real float
 * indexed
 */
#define GOVERNOR_ROPEN "GOVERNOR_ROPEN"

/**
 * Governor RCLOSE
 * type: real float
 * indexed
 */
#define GOVERNOR_RCLOSE "GOVERNOR_RCLOSE"

/**
 * Governor KIMW
 * type: real float
 * indexed
 */
#define GOVERNOR_KIMW "GOVERNOR_KIMW"

/**
 * Governor ASET
 * type: real float
 * indexed
 */
#define GOVERNOR_ASET "GOVERNOR_ASET"

/**
 * Governor KA
 * type: real float
 * indexed
 */
#define GOVERNOR_KA "GOVERNOR_KA"

/**
 * Governor TA
 * type: real float
 * indexed
 */
#define GOVERNOR_TA "GOVERNOR_TA"

/**
 * Governor TRATE
 * type: real float
 * indexed
 */
#define GOVERNOR_TRATE "GOVERNOR_TRATE"

/**
 * Governor DB
 * type: real float
 * indexed
 */
#define GOVERNOR_DB "GOVERNOR_DB"

/**
 * Governor TSA
 * type: real float
 * indexed
 */
#define GOVERNOR_TSA "GOVERNOR_TSA"

/**
 * Governor TSB
 * type: real float
 * indexed
 */
#define GOVERNOR_TSB "GOVERNOR_TSB"

/**
 * Governor RUP
 * type: real float
 * indexed
 */
#define GOVERNOR_RUP "GOVERNOR_RUP"

/**
 * Governor RDOWN
 * type: real float
 * indexed
 */
#define GOVERNOR_RDOWN "GOVERNOR_RDOWN"

/**
 * Governor TD
 * type: real float
 * indexed
 */
#define GOVERNOR_TD "GOVERNOR_TD"

/**
 * Governor KI
 * type: real float
 * indexed
 */
#define GOVERNOR_KI "GOVERNOR_KI"

/**
 * Governor TF
 * type: real float
 * indexed
 */
#define GOVERNOR_TF "GOVERNOR_TF"

/**
 * Governor KD
 * type: real float
 * indexed
 */
#define GOVERNOR_KD "GOVERNOR_KD"

/**
 * Governor KP
 * type: real float
 * indexed
 */
#define GOVERNOR_KP "GOVERNOR_KP"

/**
 * Governor TT
 * type: real float
 * indexed
 */
#define GOVERNOR_TT "GOVERNOR_TT"

/**
 * Governor GV5
 * type: real float
 * indexed
 */
#define GOVERNOR_KG "GOVERNOR_KG"

/**
 * Governor TP
 * type: real float
 * indexed
 */
#define GOVERNOR_TP "GOVERNOR_TP"

/**
 * Governor VELOPEN
 * type: real float
 * indexed
 */
#define GOVERNOR_VELOPEN "GOVERNOR_VELOPEN"

/**
 * Governor VELCLOSE
 * type: real float
 * indexed
 */
#define GOVERNOR_VELCLOSE "GOVERNOR_VELCLOSE"

/**
 * Governor ATURB
 * type: real float
 * indexed
 */
#define GOVERNOR_ATURB "GOVERNOR_ATURB"

/**
 * Governor BTURB
 * type: real float
 * indexed
 */
#define GOVERNOR_BTURB "GOVERNOR_BTURB"

/**
 * Governor TTURB
 * type: real float
 * indexed
 */
#define GOVERNOR_TTURB "GOVERNOR_TTURB"

/**
 * Governor TRATE
 * type: real float
 * indexed
 */
#define GOVERNOR_TRATE "GOVERNOR_TRATE"

// EXCITER DATA
/**
 * Flag to indicate that exciter is present
 * type: boolean
 * indexed
 */
#define HAS_EXCITER "HAS_EXCITER"

/**
 * Exciter model
 * type: string
 * indexed
 */
#define EXCITER_MODEL "EXCITER_MODEL"

/**
 * Exciter TR
 * type: real float
 * indexed
 */
#define EXCITER_TR "EXCITER_TR"

/**
 * Exciter KA
 * type: real float
 * indexed
 */
#define EXCITER_KA "EXCITER_KA"

/**
 * Exciter TA
 * type: real float
 * indexed
 */
#define EXCITER_TA "EXCITER_TA"

/**
 * Exciter TB
 * type: real float
 * indexed
 */
#define EXCITER_TB "EXCITER_TB"

/**
 * Exciter TC
 * type: real float
 * indexed
 */
#define EXCITER_TC "EXCITER_TC"

/**
 * Exciter VRMAX
 * type: real float
 * indexed
 */
#define EXCITER_VRMAX "EXCITER_VRMAX"

/**
 * Exciter VRMIN
 * type: real float
 * indexed
 */
#define EXCITER_VRMIN "EXCITER_VRMIN"

/**
 * Exciter KE
 * type: real float
 * indexed
 */
#define EXCITER_KE "EXCITER_KE"

/**
 * Exciter TE
 * type: real float
 * indexed
 */
#define EXCITER_TE "EXCITER_TE"

/**
 * Exciter KF
 * type: real float
 * indexed
 */
#define EXCITER_KF "EXCITER_KF"

/**
 * Exciter TF1
 * type: real float
 * indexed
 */
#define EXCITER_TF1 "EXCITER_TF1"

/**
 * Exciter SWITCH
 * type: real float
 * indexed
 */
#define EXCITER_SWITCH "EXCITER_SWITCH"

/**
 * Exciter E1
 * type: real float
 * indexed
 */
#define EXCITER_E1 "EXCITER_E1"

/**
 * Exciter SE1
 * type: real float
 * indexed
 */
#define EXCITER_SE1 "EXCITER_SE1"

/**
 * Exciter E2
 * type: real float
 * indexed
 */
#define EXCITER_E2 "EXCITER_E2"

/**
 * Exciter SE2
 * type: real float
 * indexed
 */
#define EXCITER_SE2 "EXCITER_SE2"

/**
 * Exciter UEL
 * type: real float
 * indexed
 */
#define EXCITER_UEL "EXCITER_UEL"

/**
 * Exciter VOS
 * type: real float
 * indexed
 */
#define EXCITER_VOS "EXCITER_VOS"

/**
 * Exciter VIMAX
 * type: real float
 * indexed
 */
#define EXCITER_VIMAX "EXCITER_VIMAX"

/**
 * Exciter VIMIN
 * type: real float
 * indexed
 */
#define EXCITER_VIMIN "EXCITER_VIMIN"

/**
 * Exciter TC1
 * type: real float
 * indexed
 */
#define EXCITER_TC1 "EXCITER_TC1"

/**
 * Exciter TB1
 * type: real float
 * indexed
 */
#define EXCITER_TB1 "EXCITER_TB1"

/**
 * Exciter VAMAX
 * type: real float
 * indexed
 */
#define EXCITER_VAMAX "EXCITER_VAMAX"

/**
 * Exciter VAMIN
 * type: real float
 * indexed
 */
#define EXCITER_VAMIN "EXCITER_VAMIN"

/**
 * Exciter KC
 * type: real float
 * indexed
 */
#define EXCITER_KC "EXCITER_KC"

/**
 * Exciter TF
 * type: real float
 * indexed
 */
#define EXCITER_TF "EXCITER_TF"

/**
 * Exciter KLR
 * type: real float
 * indexed
 */
#define EXCITER_KLR "EXCITER_KLR"

/**
 * Exciter ILR
 * type: real float
 * indexed
 */
#define EXCITER_ILR "EXCITER_ILR"

/**
 * Exciter KPR
 * type: real float
 * indexed
 */
#define EXCITER_KPR "EXCITER_KPR"

/**
 * Exciter KIR
 * type: real float
 * indexed
 */
#define EXCITER_KIR "EXCITER_KIR"

/**
 * Exciter KPM
 * type: real float
 * indexed
 */
#define EXCITER_KPM "EXCITER_KPM"

/**
 * Exciter KIM
 * type: real float
 * indexed
 */
#define EXCITER_KIM "EXCITER_KIM"

/**
 * Exciter VMMAX
 * type: real float
 * indexed
 */
#define EXCITER_VMMAX "EXCITER_VMMAX"

/**
 * Exciter VMMIN
 * type: real float
 * indexed
 */
#define EXCITER_VMMIN "EXCITER_VMMIN"

/**
 * Exciter KG
 * type: real float
 * indexed
 */
#define EXCITER_KG "EXCITER_KG"

/**
 * Exciter KP
 * type: real float
 * indexed
 */
#define EXCITER_KP "EXCITER_KP"

/**
 * Exciter KI
 * type: real float
 * indexed
 */
#define EXCITER_KI "EXCITER_KI"

/**
 * Exciter VBMAX
 * type: real float
 * indexed
 */
#define EXCITER_VBMAX "EXCITER_VBMAX"

/**
 * Exciter KC
 * type: real float
 * indexed
 */
#define EXCITER_KC "EXCITER_KC"

/**
 * Exciter XL
 * type: real float
 * indexed
 */
#define EXCITER_XL "EXCITER_XL"

/**
 * Exciter THETAP
 * type: real float
 * indexed
 */
#define EXCITER_THETAP "EXCITER_THETAP"

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
 * Alpha-numeric identifier assigned to new branches that may be
 * as part of the parsing process
 * type: string
 */
#define NEW_BRANCH_TYPE "NEW_BRANCH_TYPE"

/**
 * Number of transmission elements on branch
 * type: integer
 */
#define BRANCH_NUM_ELEMENTS "BRANCH_NUM_ELEMENTS"

/**
 * Logical flag that defines whether or not to and from bus are switched for
 * transmission element
 * type: boolean
 * indexed
 */
#define BRANCH_SWITCHED "BRANCH_SWITCHED"

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
 * Flag that indicates that branch was generated from a 3-winding transformer
 * type: boolean
 * indexed
 */
#define BRANCH_3WINDING "BRANCH_3WINDING"

/**
 * Initial branch status
 * 1: in-service
 * 0: out-of-service
 * type: integer
 * indexed
 */
#define BRANCH_STATUS "BRANCH_STATUS"

/**
 * Flag to indicate what end of the branch is metered
 * <=1: From bus is metered
 * >=2: To bus is metered
 * type: integer
 * indexed
 */
#define BRANCH_METER "BRANCH_METER"

/**
 * Parameter describing length of line
 * type: real float
 * indexed
 */
#define BRANCH_LENGTH "BRANCH_LENGTH"

/**
 * Branch owner number 1 (out of up to four)
 * type: integer
 * indexed
 */
#define BRANCH_O1 "BRANCH_O1"

/**
 * Fraction owned by owner 1
 * type: real float
 * indexed
 */
#define BRANCH_F1 "BRANCH_F1"

/**
 * Branch owner number 2 (out of up to four)
 * type: integer
 * indexed
 */
#define BRANCH_O2 "BRANCH_O2"

/**
 * Fraction owned by owner 2
 * type: real float
 * indexed
 */
#define BRANCH_F2 "BRANCH_F2"

/**
 * Branch owner number 3 (out of up to four)
 * type: integer
 * indexed
 */
#define BRANCH_O3 "BRANCH_O3"

/**
 * Fraction owned by owner 3
 * type: real float
 * indexed
 */
#define BRANCH_F3 "BRANCH_F3"

/**
 * Branch owner number 4 (out of up to four)
 * type: integer
 * indexed
 */
#define BRANCH_O4 "BRANCH_O4"

/**
 * Fraction owned by owner 4
 * type: real float
 * indexed
 */
#define BRANCH_F4 "BRANCH_F4"


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

// AREA DATA
/**
 * Total number of area fields
 * type: integer
 */
#define AREA_TOTAL "AREA_TOTAL"

/**
 * Area number
 * type: integer
 * indexed
 */
#define AREAINTG_NUMBER "AREAINTG_NUMBER"

/**
 * Bus number of the area slack bus
 * type: integer
 * indexed
 */
#define AREAINTG_ISW "AREAINTG_ISW"

/**
 * Desired net interchange leaving the area (export); entered in MW
 * type: real float
 * indexed
 */
#define AREAINTG_PDES "AREAINTG_PDES"

/**
 * Interchange tolerance bandwidth; entered in MW
 * type: real float
 * indexed
 */
#define AREAINTG_PTOL "AREAINTG_PTOL"

/**
 * Area name
 * type: string
 * indexed
 */
#define AREAINTG_NAME "AREAINTG_NAME"


// SWITCHED DATA

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
#define XFMR_CORR_TABLE_T1  "XFMR_CORR_TABLE_T1"
#define XFMR_CORR_TABLE_T2  "XFMR_CORR_TABLE_T2"
#define XFMR_CORR_TABLE_T3  "XFMR_CORR_TABLE_T3"
#define XFMR_CORR_TABLE_T4  "XFMR_CORR_TABLE_T4"
#define XFMR_CORR_TABLE_T5  "XFMR_CORR_TABLE_T5"
#define XFMR_CORR_TABLE_T6  "XFMR_CORR_TABLE_T6"
#define XFMR_CORR_TABLE_T7  "XFMR_CORR_TABLE_T7"
#define XFMR_CORR_TABLE_T8  "XFMR_CORR_TABLE_T8"
#define XFMR_CORR_TABLE_T9  "XFMR_CORR_TABLE_T9"
#define XFMR_CORR_TABLE_T10 "XFMR_CORR_TABLE_T10"
#define XFMR_CORR_TABLE_T11 "XFMR_CORR_TABLE_T11"

/**
 * Scaling factor by which transformer nominal impedance is to be multiplied to obtain
 * the actual transformer impedance for the corresponding "Ti"
 * type: real float
 */
#define XFMR_CORR_TABLE_F1  "XFMR_CORR_TABLE_F1"
#define XFMR_CORR_TABLE_F2  "XFMR_CORR_TABLE_F2"
#define XFMR_CORR_TABLE_F3  "XFMR_CORR_TABLE_F3"
#define XFMR_CORR_TABLE_F4  "XFMR_CORR_TABLE_F4"
#define XFMR_CORR_TABLE_F5  "XFMR_CORR_TABLE_F5"
#define XFMR_CORR_TABLE_F6  "XFMR_CORR_TABLE_F6"
#define XFMR_CORR_TABLE_F7  "XFMR_CORR_TABLE_F7"
#define XFMR_CORR_TABLE_F8  "XFMR_CORR_TABLE_F8"
#define XFMR_CORR_TABLE_F9  "XFMR_CORR_TABLE_F9"
#define XFMR_CORR_TABLE_F10 "XFMR_CORR_TABLE_F10"
#define XFMR_CORR_TABLE_T11 "XFMR_CORR_TABLE_T11"


// MULTISECTION LINE GROUPING 
/**
 * Two-character upper case alphanumeric multisection line grouping identifier.
 * type: string
 */
#define MULTI_SEC_LINE_ID "MULTI_SEC_LINE_ID"

/**
 * Metered end flag
 * <= 1 to designate from bus as the metered end
 * >= 2 to designate to bus as metered end
 * type: integer
 */
#define MULTI_SEC_LINE_MET "MULTI_SEC_LINE_MET"

/**
 * Bus numbers of the "dummy buses" connected by the branches that comprise this multisection
 * line grouping. No defaults allowed.
 * type: integer
 */
#define MULTI_SEC_LINE_DUM1 "MULTI_SEC_LINE_DUM1"
#define MULTI_SEC_LINE_DUM2 "MULTI_SEC_LINE_DUM2"
#define MULTI_SEC_LINE_DUM3 "MULTI_SEC_LINE_DUM3"
#define MULTI_SEC_LINE_DUM4 "MULTI_SEC_LINE_DUM4"
#define MULTI_SEC_LINE_DUM5 "MULTI_SEC_LINE_DUM5"
#define MULTI_SEC_LINE_DUM6 "MULTI_SEC_LINE_DUM6"
#define MULTI_SEC_LINE_DUM7 "MULTI_SEC_LINE_DUM7"
#define MULTI_SEC_LINE_DUM8 "MULTI_SEC_LINE_DUM8"
#define MULTI_SEC_LINE_DUM9 "MULTI_SEC_LINE_DUM9"


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

#endif /* DICTIONARY_HPP_ */
