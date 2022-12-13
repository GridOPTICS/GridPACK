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

#include "variable_defs/bus_defs.hpp"
#include "variable_defs/shunt_defs.hpp"
#include "variable_defs/load_defs.hpp"
#include "variable_defs/generator_defs.hpp"
#include "variable_defs/governor_defs.hpp"
#include "variable_defs/exciter_defs.hpp"
#include "variable_defs/branch_defs.hpp"
#include "variable_defs/transformer_defs.hpp"
#include "variable_defs/relay_defs.hpp"
#include "variable_defs/psssim_defs.hpp"
#include "variable_defs/wind_defs.hpp"
#include "variable_defs/misc_defs.hpp"

#endif /* DICTIONARY_HPP_ */
