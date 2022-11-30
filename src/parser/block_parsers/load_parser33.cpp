/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * load_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "load_parser33.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to
 internal indices
 */
gridpack::parser::LoadParser33::LoadParser33(
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
gridpack::parser::LoadParser33::~LoadParser33(void)
{
}

/**
 * parse load block
 * @param stream input stream that feeds lines from RAW file
 * @param data vector of data collection objects to store parameters
 *             from RAW file for loads
 */
void gridpack::parser::LoadParser33::parse(
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

    // LOAD_BUSNUMBER               "I"                   integer
    int l_idx, o_idx;
    o_idx = getBusIndex(split_line[0]);
    if (o_idx < 0) {
      stream.nextLine(line);
      continue;
    }
    std::map<int, int>::iterator it;
    it = p_busMap.find(o_idx);
    if (it != p_busMap.end()) {
      l_idx = it->second;
    } else {
      stream.nextLine(line);
      continue;
    }
    int nstr = split_line.size();
    // Find out how many loads are already on bus
    int nld;
    if (!p_busData[l_idx]->getValue(LOAD_NUMBER, &nld)) nld = 0;


    p_busData[l_idx]->addValue(LOAD_BUSNUMBER, o_idx, nld);

    gridpack::utility::StringUtils util;
    if (nstr > 1) {
      // Clean up 2 character tag
      std::string tag = util.clean2Char(split_line[1]);
      // LOAD_ID              "ID"                  integer
      p_busData[l_idx]->addValue(LOAD_ID, tag.c_str(), nld);
    }

    // LOAD_STATUS              "ID"                  integer
    if (nstr > 2) p_busData[l_idx]->addValue(LOAD_STATUS,
        atoi(split_line[2].c_str()), nld);

    // LOAD_AREA            "AREA"                integer
    if (nstr > 3) p_busData[l_idx]->addValue(LOAD_AREA,
        atoi(split_line[3].c_str()), nld);

    // LOAD_ZONE            "ZONE"                integer
    if (nstr > 4) p_busData[l_idx]->addValue(LOAD_ZONE,
        atoi(split_line[4].c_str()), nld);

    // LOAD_PL              "PL"                  float
    if (nstr > 5) {
      if (nld == 0) p_busData[l_idx]->addValue(LOAD_PL, atof(split_line[5].c_str()));
      p_busData[l_idx]->addValue(LOAD_PL, atof(split_line[5].c_str()), nld);
    }

    // LOAD_QL              "QL"                  float
    if (nstr > 6) {
      if (nld == 0) p_busData[l_idx]->addValue(LOAD_QL, atof(split_line[6].c_str()));
      p_busData[l_idx]->addValue(LOAD_QL, atof(split_line[6].c_str()), nld);
    }

    // LOAD_IP              "IP"                  float
    if (nstr > 7) p_busData[l_idx]->addValue(LOAD_IP,
        atof(split_line[7].c_str()), nld);

    // LOAD_IQ              "IQ"                  float
    if (nstr > 8) p_busData[l_idx]->addValue(LOAD_IQ,
        atof(split_line[8].c_str()), nld);

    // LOAD_YP              "YP"                  float
    if (nstr > 9) p_busData[l_idx]->addValue(LOAD_YP,
        atof(split_line[9].c_str()), nld);

    // LOAD_YQ            "YQ"                integer
    if (nstr > 10) p_busData[l_idx]->addValue(LOAD_YQ,
        atof(split_line[10].c_str()), nld);

    // TODO: add variables OWNER, SCALE, INTRPT

    // Increment number of loads in data object
    if (nld == 0) {
      nld = 1;
      p_busData[l_idx]->addValue(LOAD_NUMBER,nld);
    } else {
      nld++;
      p_busData[l_idx]->setValue(LOAD_NUMBER,nld);
    }

    stream.nextLine(line);
  }
}
