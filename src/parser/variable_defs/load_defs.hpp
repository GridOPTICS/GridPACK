/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all load parameters that can be read in
 * from PTI format files. Each parameter has a corresponding macro that can be
 * used as a unique string to identify the parameter. The use of macros instead
 * of using strings directly will provide extra safety by forcing compiler errors
 * in the case of typos or spelling mistakes.
 */

/**
 *  Variables that can be associated more than once for a bus or a branch can be
 *  indexed by an integer to distinguish different instances. For example,
 *  multiple generators can be associated with a bus and multiple transmission
 *  elements can be associated with a branch. The variables that have an associated
 *  index are denoted with the keyword "indexed".
 */

#ifndef _LOAD_VAR_HPP_
#define _LOAD_VAR_HPP_

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
 * Load DYN_PERC
 * type: float
 * indexed
 */
#define LOAD_DYN_PERC "LOAD_DYN_PERC"

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
 * LOAD_AC_PERC, ac motor load percentage
 * type: float
 * indexed
 */
#define LOAD_AC_PERC "LOAD_AC_PERC"

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

#endif /* _LOAD_VAR_HPP_ */
