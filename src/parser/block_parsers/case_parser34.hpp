/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * case_parser34.hpp
 *       Created on: December 8, 2022
 *           Author: Bruce Palmer
 */
#ifndef _CASE_PARSER34_H
#define _CASE_PARSER34_H

#include "gridpack/parser/block_parsers/base_block_parser.hpp"

namespace gridpack {
namespace parser {

class CaseParser34 : public BaseBlockParser {
  public:
  /**
   * Constructor
   * @param bus_map map indices in RAW file to internal indices
   * @param name_map map name in RAW file to internal indices
   * @param branch_map map bus index pair in RAW file to internal indices
   */
  CaseParser34(
      std::map<int,int> *bus_map,
      std::map<std::string,int> *name_map,
      std::map<std::pair<int, int>, int> *branch_map);

  /**
   * Simple Destructor
   */
  virtual ~CaseParser34(void);

  /**
   * parse case block
   * @param stream input stream that feeds lines from RAW file
   * @param data data collection object to store parameters from RAW file
   * @param sbase value of sbase from RAW file
   * @param id value if id from RAW file
   */
  void parse(
      gridpack::stream::InputStream &stream,
      boost::shared_ptr<gridpack::component::DataCollection> &data,
      double &sbase, int &id);
};

} // parser
} // gridpack
#endif
