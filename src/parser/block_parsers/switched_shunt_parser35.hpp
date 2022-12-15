/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * switched_shunt_parser35.hpp
 *       Created on: December 6, 2022
 *           Author: Bruce Palmer
 */
#ifndef _SWITCHED_SHUNT_PARSER35_H
#define _SWITCHED_SHUNT_PARSER35_H

#include "gridpack/parser/block_parsers/base_block_parser.hpp"

namespace gridpack {
namespace parser {

class SwitchedShuntParser35 : public BaseBlockParser {
  public:
  /**
   * Constructor
   * @param bus_map map indices in RAW file to internal indices
   * @param name_map map name in RAW file to internal indices
   * @param branch_map map bus index pair in RAW file to internal indices
   */
  SwitchedShuntParser35(
      std::map<int,int> *bus_map,
      std::map<std::string,int> *name_map,
      std::map<std::pair<int, int>, int> *branch_map);

  /**
   * Simple Destructor
   */
  virtual ~SwitchedShuntParser35(void);

  /**
   * parse switched shunt block
   * @param stream input stream that feeds lines from RAW file
   * @param data vector of data collection objects to store parameters
   *             from RAW file for switched shunts
   */
  void parse(
      gridpack::stream::InputStream &stream,
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &data);
};

} // parser
} // gridpack
#endif
