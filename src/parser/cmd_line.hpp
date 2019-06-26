// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   cmd_line.hpp
 * @author Bruce Palmer
 * @date   2019-05-24 09:25:03 d3g293
 * 
 * @brief  
 * This utility is designed to support parsing of command line arguments coming
 * in via the standard (int argv, char **argc) arguments of main. This module
 * ignores the first argument (the program name). The number of arguments
 * returned by the class is argv-1.
 * 
 */

// -------------------------------------------------------------

#ifndef _cmd_line_hpp_
#define _cmd_line_hpp_

#include <string>
#include <vector>
#include "gridpack/utilities/string_utils.hpp"

namespace gridpack {
namespace parser {

// -------------------------------------------------------------
//  class CmdLineParser
// -------------------------------------------------------------
class CmdLineParser {
public:

  /**
   * Default constructor
   * @param argv standard argv from main
   * @param argc standard argc from main
   */
  CmdLineParser(int argv, char** argc);

  /**
   * Default destructor
   */
  ~CmdLineParser(void);

  /**
   * Return the total number of tokens found in command line arguments. This is
   * one less than the original argv
   * @return number of command line tokens
   */
  int numCmdTokens();

  /**
   * Return the requested token (by index)
   * @param idx index of requested token
   * @return requested token
   */
  std::string getToken(int idx);

  /**
   * Get the token after the token represented by the string. The token given in
   * the argument will match any token that contains the same first characters
   * as the input string
   * @param key string that is used to match token in internal list
   * @return next token in command list after key
   */
  std::string getTokenAfter(const char *key);

private:

  // Number of tokens in command line arguments
  int p_argv;

  // Individual tokens in command line arguments
  std::vector<std::string> p_argc;
};


} // namespace gridpack
} // namespace parser

#endif
