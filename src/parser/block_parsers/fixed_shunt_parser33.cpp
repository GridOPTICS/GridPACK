/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * fixed_shunt_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "fixed_shunt_parser33.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to
 internal indices
 */
gridpack::parser::FixedShuntParser33::FixedShuntParser33(
    std::map<int,int> &bus_map,
    std::map<std::string,int> &name_map,
    std::map<std::pair<int, int>, int> &branch_map) :
    gridpack::parser::BaseBlockParser(
      bus_map, name_map, branch_map)
{
  p_busMap = bus_map;
  p_nameMap = name_map;
  p_branchMap = branch_map;
}


/**
 * Simple Destructor
 */
gridpack::parser::FixedShuntParser33::~FixedShuntParser33(void)
{
}

/**
 * parse fixed shunt block
 * @param stream input stream that feeds lines from RAW file
 * @param data vector of data collection objects to store parameters
 *             from RAW file for fixed shunts
 */
void gridpack::parser::FixedShuntParser33::parse(
    gridpack::stream::InputStream &stream,
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &p_busData)
{
  std::string          line;
  stream.nextLine(line); //this should be the first line of the block

  while(test_end(line)) {
    std::vector<std::string>  split_line;
    if (check_comment(line)) {
      stream.nextLine(line);
      continue;
    }
    this->cleanComment(line);
    split_line = this->splitPSSELine(line);

    // SHUNT_BUSNUMBER               "I"                   integer
    int l_idx, o_idx;
    o_idx = getBusIndex(split_line[0]);
    std::map<int, int>::iterator it;
    it = p_busMap.find(o_idx);
    if (it != p_busMap.end()) {
      l_idx = it->second;
    } else {
      stream.nextLine(line);
      continue;
    }
    int nstr = split_line.size();

    // Find out how many shunts are already on bus
    int nshnt;
    if (!p_busData[l_idx]->getValue(SHUNT_NUMBER, &nshnt)) nshnt = 0;

    p_busData[l_idx]->addValue(SHUNT_BUSNUMBER, o_idx, nshnt);

    if (nstr > 1) {
      // Clean up 2 character tag
      gridpack::utility::StringUtils util;
      std::string tag = util.clean2Char(split_line[1]);
      // SHUNT_ID              "ID"                  integer
      p_busData[l_idx]->addValue(SHUNT_ID, tag.c_str(), nshnt);
    }

    // SHUNT_STATUS              "STATUS"                  integer
    if (nstr > 2) p_busData[l_idx]->addValue(SHUNT_STATUS,
        atoi(split_line[2].c_str()), nshnt);

    // BUS_SHUNT_GL              "GL"                  float
    if (nstr > 3) {
      if (nshnt==0) p_busData[l_idx]->addValue(BUS_SHUNT_GL,
          atof(split_line[3].c_str()));
      p_busData[l_idx]->addValue(BUS_SHUNT_GL,
          atof(split_line[3].c_str()),nshnt);
    }

    // BUS_SHUNT_BL              "BL"                  float
    if (nstr > 4) {
      if (nshnt == 0) p_busData[l_idx]->addValue(BUS_SHUNT_BL,
          atof(split_line[4].c_str()));
      p_busData[l_idx]->addValue(BUS_SHUNT_BL,
          atof(split_line[4].c_str()),nshnt);
    }

    // Increment number of shunts in data object
    if (nshnt == 0) {
      nshnt = 1;
      p_busData[l_idx]->addValue(SHUNT_NUMBER,nshnt);
    } else {
      nshnt++;
      p_busData[l_idx]->setValue(SHUNT_NUMBER,nshnt);
    }

    stream.nextLine(line);
  }
}
