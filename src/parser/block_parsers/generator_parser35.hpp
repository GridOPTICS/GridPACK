/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * generator_parser35.hpp
 *       Created on: December 5, 2022
 *           Author: Bruce Palmer
 */
#ifndef _GENERATOR_PARSER35_H
#define _GENERATOR_PARSER35_H

#include "gridpack/parser/block_parsers/base_block_parser.hpp"

namespace gridpack {
namespace parser {

class GeneratorParser35 : public BaseBlockParser {
  public:
  /**
   * Constructor
   * @param bus_map map indices in RAW file to internal indices
   * @param name_map map name in RAW file to internal indices
   * @param branch_map map bus index pair in RAW file to internal indices
   */
  GeneratorParser35(
      std::map<int,int> *bus_map,
      std::map<std::string,int> *name_map,
      std::map<std::pair<int, int>, int> *branch_map);

  /**
   * Simple Destructor
   */
  virtual ~GeneratorParser35(void);

  /**
   * parse generator block
   * @param stream input stream that feeds lines from RAW file
   * @param data vector of data collection objects to store parameters
   *             from RAW file for generators
   */
  void parse(
      gridpack::stream::InputStream &stream,
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &data);
};

} // parser
} // gridpack
#endif
