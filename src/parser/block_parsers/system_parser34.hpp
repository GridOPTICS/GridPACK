/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * system_parser33.hpp
 *       Created on: December 5, 2022
 *           Author: Bruce Palmer
 */
#ifndef _SYSTEM_PARSER35_H
#define _SYSTEM_PARSER35_H

#include "gridpack/parser/block_parsers/base_block_parser.hpp"

namespace gridpack {
namespace parser {

class SystemParser35 : public BaseBlockParser {
  public:
  /**
   * Constructor
   * @param bus_map map indices in RAW file to internal indices
   * @param name_map map name in RAW file to internal indices
   * @param branch_map map bus index pair in RAW file to internal indices
   */
  SystemParser35(
      std::map<int,int> *bus_map,
      std::map<std::string,int> *name_map,
      std::map<std::pair<int, int>, int> *branch_map);

  /**
   * Simple Destructor
   */
  virtual ~SystemParser35(void);

  /**
   * parse System block. Currently does not store data
   * @param stream input stream that feeds lines from RAW file
   */
  void parse(
      gridpack::stream::InputStream &stream);
};

} // parser
} // gridpack
#endif
