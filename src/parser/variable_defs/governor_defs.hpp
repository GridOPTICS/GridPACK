/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all governor parameters that can be read
 * in from PTI format files. Each parameter has a corresponding macro that can
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

#ifndef _GOVERNOR_VAR_HPP_
#define _GOVERNOR_VAR_HPP_

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
 * Governor TR
 * type: real float
 * indexed
 */
#define GOVERNOR_TR "GOVERNOR_TR"

/**
 * Governor TG
 * type: real float
 * indexed
 */
#define GOVERNOR_TG "GOVERNOR_TG"

/**
 * Governor TW
 * type: real float
 * indexed
 */
#define GOVERNOR_TW "GOVERNOR_TW"


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
 * Governor r
 * type: real float
 * indexed
 */
#define GOVERNOR_r "GOVERNOR_r"

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
 * Governor VELM
 * type: real float
 * indexed
 */
#define GOVERNOR_VELM "GOVERNOR_VELM"

/**
 * Governor GMAX
 * type: real float
 * indexed
 */
#define GOVERNOR_GMAX "GOVERNOR_GMAX"

/**
 * Governor GMIN
 * type: real float
 * indexed
 */
#define GOVERNOR_GMIN "GOVERNOR_GMIN"

/**
 * Governor QNL
 * type: real float
 * indexed
 */
#define GOVERNOR_QNL "GOVERNOR_QNL"


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
 * Governor DT
 * type: real float
 * indexed
 */
#define GOVERNOR_DT "GOVERNOR_DT"

/**
 * Governor AT
 * type: real float
 * indexed
 */
#define GOVERNOR_AT "GOVERNOR_AT"
/**
 * Governor KT
 * type: real float
 * indexed
 */
#define GOVERNOR_KT "GOVERNOR_KT"


/**
 * Governor DT
 * type: real float
 * indexed
 */
#define GOVERNOR_DT "GOVERNOR_DT"

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

#endif /* DICTIONARY_HPP_ */
