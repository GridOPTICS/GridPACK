/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * imped_corr_parser33.hpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#ifndef _IMPED_CORR_PARSER33_H
#define _IMPED_CORR_PARSER33_H

#include "gridpack/parser/block_parsers/base_block_parser.hpp"

namespace gridpack {
namespace parser {

class ImpedCorrParser33 : public BaseBlockParser {
  public:
  /**
   * Constructor
   * @param bus_map map indices in RAW file to internal indices
   * @param name_map map name in RAW file to internal indices
   * @param branch_map map bus index pair in RAW file to internal indices
   */
  ImpedCorrParser33(
      std::map<int,int> *bus_map,
      std::map<std::string,int> *name_map,
      std::map<std::pair<int, int>, int> *branch_map);

  /**
   * Simple Destructor
   */
  virtual ~ImpedCorrParser33(void);

  /**
   * parse impedence correction block
   * @param stream input stream that feeds lines from RAW file
   * @param data vector of data collection objects to store parameters
   *             from RAW file for impedence correction table
   */
  void parse(
      gridpack::stream::InputStream &stream,
      std::map<int,boost::shared_ptr<gridpack::component::DataCollection> > &data);

};

} // parser
} // gridpack
#endif
