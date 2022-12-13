/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * This file contains definitions for all branch parameters that can be read in
 * from PTI format files. Each parameter has a corresponding macro that can be used
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

#ifndef _BRANCH_VAR_HPP_
#define _BRANCH_VAR_HPP_

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
 * Alphanumeric string assigned to branch
 * type: string
 * indexed
 */
#define BRANCH_NAME "BRANCH_NAME"

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
 * Nth rating for branch
 * type: real float
 * indexed
 */
#define BRANCH_RATE1 "BRANCH_RATE1"
#define BRANCH_RATE2 "BRANCH_RATE2"
#define BRANCH_RATE3 "BRANCH_RATE3"
#define BRANCH_RATE4 "BRANCH_RATE4"
#define BRANCH_RATE5 "BRANCH_RATE5"
#define BRANCH_RATE6 "BRANCH_RATE6"
#define BRANCH_RATE7 "BRANCH_RATE7"
#define BRANCH_RATE8 "BRANCH_RATE8"
#define BRANCH_RATE9 "BRANCH_RATE9"
#define BRANCH_RATE10 "BRANCH_RATE10"
#define BRANCH_RATE11 "BRANCH_RATE11"
#define BRANCH_RATE12 "BRANCH_RATE12"

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

// Branch sequence data
/**
 * Zero sequence branch resistance
 * type: real float
 * indexed
 */
#define BRANCH_SEQ_RLINZ "BRANCH_SEQ_RLINZ"

/**
 * Zero sequence branch reactance
 * type: real float
 * indexed
 */
#define BRANCH_SEQ_XLINZ "BRANCH_SEQ_XLINZ"

/**
 * Total zero sequence branch branch charging susceptance
 * type: real float
 * indexed
 */
#define BRANCH_SEQ_BCHZ "BRANCH_SEQ_BCHZ"

/**
 * Real zero sequence admittance of the line connect to "from" bus
 * type: real float
 * indexed
 */
#define BRANCH_SEQ_GI "BRANCH_SEQ_GI"

/**
 * Imaginary zero sequence admittance of the line connect to "from" bus
 * type: real float
 * indexed
 */
#define BRANCH_SEQ_BJ "BRANCH_SEQ_BJ"

/**
 * Real zero sequence admittance of the line connect to "to" bus
 * type: real float
 * indexed
 */
#define BRANCH_SEQ_GJ "BRANCH_SEQ_GJ"

/**
 * Imaginary zero sequence admittance of the line connect to "to" bus
 * type: real float
 * indexed
 */
#define BRANCH_SEQ_BI "BRANCH_SEQ_BI"

#endif /* _BRANCH_VAR_HPP_ */
