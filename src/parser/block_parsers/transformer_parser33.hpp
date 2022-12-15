/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * transformer_parser33.hpp
 *       Created on: December 1, 2022
 *           Author: Bruce Palmer
 */
#ifndef _TRANSFORMER_PARSER33_H
#define _TRANSFORMER_PARSER33_H

#include "gridpack/parser/block_parsers/base_block_parser.hpp"

namespace gridpack {
namespace parser {

class TransformerParser33 : public BaseBlockParser {
  public:
  /**
   * Constructor
   * @param bus_map map indices in RAW file to internal indices
   * @param name_map map name in RAW file to internal indices
   * @param branch_map map bus index pair in RAW file to internal indices
   */
  TransformerParser33(
      std::map<int,int> *bus_map,
      std::map<std::string,int> *name_map,
      std::map<std::pair<int, int>, int> *branch_map);

  /**
   * Simple Destructor
   */
  virtual ~TransformerParser33(void);

  /**
   * parse transformer block
   * @param stream input stream that feeds lines from RAW file
   * @param busData vector of bus data collection objects to store parameters
   *             from RAW file for transformers
   * @param branchData vector of branch data collection objects to store
   *             parameters from RAW file for transformers
   * @param case_sbase sbase parameter
   * @param maxBusIndex maximum value of bus index for any bus in system
   */
  void parse(
      gridpack::stream::InputStream &stream,
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &busData,
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &branchData,
      double case_sbase,
      int maxBusIndex);


  private:

  /**
   * Utility function to parse lines in 3-winding transformer block
   * @param split_line list of tokens from parsing a single line
   * @param windv winding parameter
   * @param ang angle parameter
   * @param ratea, rateb, ratec rating parameters
   * @return true if all parameters found
   */
  bool parse3WindXForm(std::vector<std::string>split_line, double *windv,
      double* ang, double *ratea, double *rateb, double *ratec)
  {
    *windv = atof(split_line[0].c_str());
    *ang = atof(split_line[2].c_str());
    *ratea = atof(split_line[3].c_str());
    *rateb = atof(split_line[4].c_str());
    *ratec = atof(split_line[5].c_str());
    bool ret = true;
    return ret && (split_line.size() > 5);
  }
};

} // parser
} // gridpack
#endif
