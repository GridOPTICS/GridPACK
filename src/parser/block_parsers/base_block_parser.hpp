/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 * This class defines a set of utility functions that are used by all block
 * parsers. Each block parser must define a 'parse' function for its particular
 * block. Since these functions vary widely across different blocks, a prototype
 * is not defined in the base class.
 *
 * base_block_parser.hpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#ifndef _BASE_BLOCK_PARSER_H
#define _BASE_BLOCK_PARSER_H
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp> // needed of is_any_of()
#include "gridpack/stream/input_stream.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/utilities/string_utils.hpp"
#include "gridpack/parser/dictionary.hpp"
#include <string>
#include <map>

#define TERM_CHAR '0'

namespace gridpack {
namespace parser {

class BaseBlockParser{
  public:
  /**
   * Constructor
   * @param bus_map map indices in RAW file to internal indices
   * @param name_map map name in RAW file to internal indices
   * @param branch_map map bus index pair in RAW file to internal indices
   */
  BaseBlockParser(
            std::map<int,int> &bus_map,
            std::map<std::string,int> &name_map,
            std::map<std::pair<int, int>, int> &branch_map);
  /**
   * Simple destructor
   */
  virtual ~BaseBlockParser();

  /**
   * Test to see if string is a comment line. Check to see if first
   * non-blank characters are "//"
   * @param str input string
   */
  bool check_comment(std::string &str) const;

  /**
   * Remove leading and trailing blanks
   * @param str input string
   */
  void remove_blanks(std::string &str);

  /**
   * Test to see if string terminates a section
   * @param str input string
   * @return: false if first non-blank character is TERM_CHAR
   */
  bool test_end(std::string &str) const;

  /**
   * Check to see if string represents a character string. Look for single
   * quotes. Strip off white space and remove quotes if returning true
   * @param str input string
   * @return true if string is a character string delimited by single quotes
   */
  bool check_string(std::string &str);

  /**
   * Split bus name into its name and base bus voltage. Assume that incoming
   * string has been stripped of single quotes and leading and trailing
   * white space. The name and the voltage are suppose to be delimited by at
   * least 3 white spaces
   * @param string input string
   * @param name bus name derived from string
   * @param voltage base bus voltage derived from string
   */
  void parseBusName(std::string &string, std::string &name, double &voltage);

  /**
   * Remove comment from string (all text after a single '/' character)
   * Check to see if '/' character occurs in a text string (delimited by
   * either '' or "")
   * @param string line of text to be cleaned
   */
  void cleanComment(std::string &string);

  /**
   * Split PSS/E formatted lines into individual tokens using both blanks and
   * commas as delimiters
   * @param line input string from PSS/E file
   * @return vector of tokens parsed from PSS/E line
   */
  std::vector<std::string> splitPSSELine (std::string line);

  /**
   * Check to see if string is blank
   * @param string string that needs to checked for non-blank characters
   * @return true if no non-blank characters are found
   */
  bool isBlank(std::string string);

private:
  std::map<int,int> p_busMap;
  std::map<std::string,int> p_nameMap;
  std::map<std::pair<int, int>, int> p_branchMap;
};

} // parse
} // gridpack
#endif
