/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * interarea_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "interarea_parser33.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to internal indices
 */
gridpack::parser::InterAreaParser33::InterAreaParser33(
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
gridpack::parser::InterAreaParser33::~InterAreaParser33(void)
{
}

/**
 * parse inter-area block. Currently does not store data
 * @param stream input stream that feeds lines from RAW file
 */
void gridpack::parser::InterAreaParser33::parse(
    gridpack::stream::InputStream &stream)
{
      std::string          line;

      stream.nextLine(line); //this should be the first line of the block

      while(test_end(line)) {
#if 0
        std::vector<std::string>  split_line;
        if (check_comment(line)) {
          stream.nextLine(line);
          continue;
        }
        split_line = this->splitPSSELine(line);
        std::vector<gridpack::component::DataCollection>   inter_area_instance;
        gridpack::component::DataCollection          data;

        /*
         * type: integer
         * #define INTERAREA_TRANSFER_FROM "INTERAREA_TRANSFER_FROM"
         */
        data.addValue(INTERAREA_TRANSFER_FROM, atoi(split_line[0].c_str()));
        inter_area_instance.push_back(data);

        /*
         * type: integer
         * #define INTERAREA_TRANSFER_TO "INTERAREA_TRANSFER_TO"
         */
        data.addValue(INTERAREA_TRANSFER_TO, atoi(split_line[0].c_str()));
        inter_area_instance.push_back(data);

        /*
         * type: character
         * #define INTERAREA_TRANSFER_TRID "INTERAREA_TRANSFER_TRID"
         */
        data.addValue(INTERAREA_TRANSFER_TRID, split_line[0].c_str()[0]);
        inter_area_instance.push_back(data);

        /*
         * type: real float
         * #define INTERAREA_TRANSFER_PTRAN "INTERAREA_TRANSFER_PTRAN"
         */
        data.addValue(INTERAREA_TRANSFER_PTRAN, atof(split_line[0].c_str()));
        inter_area_instance.push_back(data);

        inter_area.push_back(inter_area_instance);
#endif
        stream.nextLine(line);
      }
}
