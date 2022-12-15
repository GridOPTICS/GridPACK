/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * system_parser34.cpp
 *       Created on: December December 5, 2022
 *           Author: Bruce Palmer
 */
#include "system_parser34.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to internal indices
 */
gridpack::parser::SystemParser34::SystemParser34(
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
gridpack::parser::SystemParser34::~SystemParser34(void)
{
}

/**
 * parse System block. Currently does not store data
 * @param stream input stream that feeds lines from RAW file
 */
void gridpack::parser::SystemParser34::parse(
    gridpack::stream::InputStream &stream)
{
  std::string          line;

  stream.nextLine(line); //this should be the first line of the block

  while(test_end(line)) {
    // TODO: parse something here
    stream.nextLine(line);
  }
    
}
