/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * owner_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "owner_parser33.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to internal indices
 */
gridpack::parser::OwnerParser33::OwnerParser33(
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
gridpack::parser::OwnerParser33::~OwnerParser33(void)
{
}

/**
 * parse owner block
 * @param stream input stream that feeds lines from RAW file
 * @param p_network_data data collection object to store parameters from RAW file
 */
void gridpack::parser::OwnerParser33::parse(
    gridpack::stream::InputStream &stream,
    boost::shared_ptr<gridpack::component::DataCollection> &p_network_data)
{
  std::string          line;
  stream.nextLine(line); //this should be the first line of the block

  int ncnt = 0;
  while(test_end(line)) {
    std::vector<std::string>  split_line;
    if (check_comment(line)) {
      stream.nextLine(line);
      continue;
    }
    this->cleanComment(line);
    split_line = this->splitPSSELine(line);

    if (split_line.size() >= 2) {
      // OWNER_NUMBER             "I"                    integer
      p_network_data->addValue(OWNER_NUMBER, atoi(split_line[0].c_str()), ncnt);

      // OWNER_NAME              "OWNAM"                 string
      p_network_data->addValue(OWNER_NAME, split_line[1].c_str(), ncnt);
      ncnt++;
    }

    stream.nextLine(line);
  }
  p_network_data->addValue(OWNER_TOTAL, ncnt);
}
