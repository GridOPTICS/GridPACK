/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * transformer_parser34.hpp
 *       Created on: December 9, 2022
 *           Author: Bruce Palmer
 */
#ifndef _TRANSFORMER_PARSER34_H
#define _TRANSFORMER_PARSER34_H

#include "gridpack/parser/block_parsers/base_block_parser.hpp"

namespace gridpack {
namespace parser {

class TransformerParser34 : public BaseBlockParser {
  public:
  /**
   * Constructor
   * @param bus_map map indices in RAW file to internal indices
   * @param name_map map name in RAW file to internal indices
   * @param branch_map map bus index pair in RAW file to internal indices
   */
  TransformerParser34(
      std::map<int,int> *bus_map,
      std::map<std::string,int> *name_map,
      std::map<std::pair<int, int>, int> *branch_map);

  /**
   * Simple Destructor
   */
  virtual ~TransformerParser34(void);

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
      double* ang, double *rate)
  {
    *windv = atof(split_line[0].c_str());
    *ang = atof(split_line[2].c_str());
    rate[0] = atof(split_line[3].c_str());
    rate[1] = atof(split_line[4].c_str());
    rate[2] = atof(split_line[5].c_str());
    rate[3] = atof(split_line[6].c_str());
    rate[4] = atof(split_line[7].c_str());
    rate[5] = atof(split_line[8].c_str());
    rate[6] = atof(split_line[9].c_str());
    rate[7] = atof(split_line[10].c_str());
    rate[8] = atof(split_line[11].c_str());
    rate[9] = atof(split_line[12].c_str());
    rate[10] = atof(split_line[13].c_str());
    rate[11] = atof(split_line[14].c_str());
    bool ret = true;
    return ret && (split_line.size() > 14);
  }
};

} // parser
} // gridpack
#endif
