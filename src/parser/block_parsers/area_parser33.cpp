/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * area_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "area_parser33.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to internal indices
 */
gridpack::parser::AreaParser33::AreaParser33(
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
gridpack::parser::AreaParser33::~AreaParser33(void)
{
}

/**
 * parse area block
 * @param stream input stream that feeds lines from RAW file
 * @param p_network_data data collection object to store parameters from RAW file
 */
void gridpack::parser::AreaParser33::parse(
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

    // AREAINTG_ISW             "I"                    integer
    p_network_data->addValue(AREAINTG_ISW, atoi(split_line[1].c_str()),ncnt);

    // AREAINTG_NUMBER             "I"                    integer
    p_network_data->addValue(AREAINTG_NUMBER, atoi(split_line[0].c_str()),ncnt);

    // AREAINTG_PDES          "PDES"                 float
    p_network_data->addValue(AREAINTG_PDES, atof(split_line[2].c_str()),ncnt);

    // AREAINTG_PTOL          "PTOL"                 float
    p_network_data->addValue(AREAINTG_PTOL, atof(split_line[3].c_str()),ncnt);

    // AREAINTG_NAME         "ARNAM"                string
    p_network_data->addValue(AREAINTG_NAME, split_line[4].c_str(),ncnt);
    ncnt++;

    stream.nextLine(line);
  }
  p_network_data->addValue(AREA_TOTAL,ncnt);
}
