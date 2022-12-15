/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all exciter parameters that can be read in
 * from PTI format files. Each parameter has a corresponding macro that can be
 * used as a unique string to identify the parameter. The use of macros instead of
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

#ifndef _EXCITER_VAR_HPP_
#define _EXCITER_VAR_HPP_

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
 * Exciter K
 * type: real float
 * indexed
 */
#define EXCITER_K "EXCITER_K"

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
 * Exciter TA_OVER_TB
 * type: real float
 * indexed
 */
#define EXCITER_TA_OVER_TB "EXCITER_TA_OVER_TB"

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
 * Exciter EMAX
 * type: real float
 * indexed
 */
#define EXCITER_EMAX "EXCITER_EMAX"

/**
 * Exciter EMIN
 * type: real float
 * indexed
 */
#define EXCITER_EMIN "EXCITER_EMIN"


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

#endif /* _EXCITER_VAR_HPP_ */
