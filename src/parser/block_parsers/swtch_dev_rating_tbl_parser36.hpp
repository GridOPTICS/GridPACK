/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 *       swtch_dev_rating_tbl_parser36.hpp
 *       Created on: February 17, 2025
 *           Author: Bruce Palmer
 */
#ifndef _SWTCH_DEV_RATING_TBL_PARSER36_H
#define _SWTCH_DEV_RATING_TBL_PARSER36_H

#include "gridpack/parser/block_parsers/base_block_parser.hpp"

namespace gridpack {
namespace parser {

class SwtchDevRatingTblParser36 : public BaseBlockParser {
  public:
  /**
   * Constructor
   * @param bus_map map indices in RAW file to internal indices
   * @param name_map map name in RAW file to internal indices
   * @param branch_map map bus index pair in RAW file to internal indices
   */
  SwtchDevRatingTblParser36(
      std::map<int,int> *bus_map,
      std::map<std::string,int> *name_map,
      std::map<std::pair<int, int>, int> *branch_map);

  /**
   * Simple Destructor
   */
  virtual ~SwtchDevRatingTblParser36(void);

  /**
   * parse interarea block. Currently does not store data
   * @param stream input stream that feeds lines from RAW file
   */
  void parse(
      gridpack::stream::InputStream &stream);
};

} // parser
} // gridpack
#endif
