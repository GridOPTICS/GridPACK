/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * vsc_line_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "vsc_line_parser33.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to internal indices
 */
gridpack::parser::VSCLineParser33::VSCLineParser33(
    std::map<int,int> *bus_map,
    std::map<std::string,int> *name_map,
    std::map<std::pair<int, int>, int> *branch_map) :
    gridpack::parser::BaseBlockParser(
      bus_map, name_map, branch_map)
{
}


/**
 * Simple Destructor
 */
gridpack::parser::VSCLineParser33::~VSCLineParser33(void)
{
}

/**
 * parse vsc line block. Currently does not store data
 * @param stream input stream that feeds lines from RAW file
 */
void gridpack::parser::VSCLineParser33::parse(
    gridpack::stream::InputStream &stream)
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
    int l_idx, o_idx;
    o_idx = atoi(split_line[1].c_str());
    std::map<int, int>::iterator it;
    it = p_busMap->find(o_idx);
    if (it != p_busMap->end()) {
      l_idx = it->second;
    } else {
      stream.nextLine(line);
      continue;
    }

    stream.nextLine(line);
  }
}
