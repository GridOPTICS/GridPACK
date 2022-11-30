/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * case_parser33.hpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#ifndef _CASE_PARSER33_H
#define _CASE_PARSER33_H

#include "gridpack/parser/block_parsers/base_block_parser.hpp"

namespace gridpack {
namespace parser {

class CaseParser33 : public BaseBlockParser {
  public:
  /**
   * Constructor
   * @param stream input stream that feeds lines from RAW file
   * @param name_map map name in RAW file to internal indices
   * @param branch_map map bus index pair in RAW file to internal indices
   */
  CaseParser33(
      std::map<int,int> &bus_map,
      std::map<std::string,int> &name_map,
      std::map<std::pair<int, int>, int> &branch_map);

  /**
   * Simple Destructor
   */
  virtual ~CaseParser33(void);

  /**
   * parse case block
   * @param bus_map map indices in RAW file to internal indices
   * @param data data collection object to store parameters from RAW file
   * @param sbase value of sbase from RAW file
   * @param id value if id from RAW file
   */
  void parse(
      gridpack::stream::InputStream &stream,
      boost::shared_ptr<gridpack::component::DataCollection> &data,
      double &sbase, int &id);

  private:

  std::map<int,int> p_busMap;
  std::map<std::string,int> p_nameMap;
  std::map<std::pair<int, int>, int> p_branchMap;

  double p_case_sbase;
  int p_case_id;

};

} // parser
} // gridpack
#endif
