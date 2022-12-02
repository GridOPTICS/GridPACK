/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * switched_shunt_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "switched_shunt_parser33.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to
 internal indices
 */
gridpack::parser::SwitchedShuntParser33::SwitchedShuntParser33(
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
gridpack::parser::SwitchedShuntParser33::~SwitchedShuntParser33(void)
{
}

/**
 * parse switched shunt block
 * @param stream input stream that feeds lines from RAW file
 * @param data vector of data collection objects to store parameters
 *             from RAW file for switched shunts
 */
void gridpack::parser::SwitchedShuntParser33::parse(
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

    /*
     * type: integer
     * #define SWSHUNT_BUSNUMBER "SWSHUNT_BUSNUMBER"
     */
    int l_idx, o_idx;
    l_idx = atoi(split_line[0].c_str());
#ifdef OLD_MAP
    std::map<int, int>::iterator it;
#else
    boost::unordered_map<int, int>::iterator it;
#endif
    it = p_busMap->find(l_idx);
    if (it != p_busMap->end()) {
      o_idx = it->second;
    } else {
      stream.nextLine(line);
      continue;
    }
    int nval = split_line.size();

    p_busData[o_idx]->addValue(SWSHUNT_BUSNUMBER, atoi(split_line[0].c_str()));

    /*
     * type: integer
     * #define SHUNT_MODSW "SHUNT_MODSW"
     */
    p_busData[o_idx]->addValue(SHUNT_MODSW, atoi(split_line[1].c_str()));

    /*
     * type: integer
     * #define SHUNT_ADJM "SHUNT_ADJM"
     */
    p_busData[o_idx]->addValue(SHUNT_ADJM, atoi(split_line[2].c_str()));

    /*
     * type: integer
     * #define SHUNT_SWCH_STAT "SHUNT_SWCH_STAT"
     */
    p_busData[o_idx]->addValue(SHUNT_SWCH_STAT, atoi(split_line[3].c_str()));

    /*
     * type: real float
     * #define SHUNT_VSWHI "SHUNT_VSWHI"
     */
    p_busData[o_idx]->addValue(SHUNT_VSWHI, atof(split_line[4].c_str()));

    /*
     * type: real float
     * #define SHUNT_VSWLO "SHUNT_VSWLO"
     */
    p_busData[o_idx]->addValue(SHUNT_VSWLO, atof(split_line[5].c_str()));

    /*
     * type: integer
     * #define SHUNT_SWREM "SHUNT_SWREM"
     */
    p_busData[o_idx]->addValue(SHUNT_SWREM, atoi(split_line[6].c_str()));

    /*
     * type: real float
     * #define SHUNT_RMPCT "SHUNT_RMPCT"
     */
    p_busData[o_idx]->addValue(SHUNT_RMPCT, atof(split_line[7].c_str()));

    /*
     * type: string
     * #define SHUNT_RMIDNT "SHUNT_RMIDNT"
     */
    p_busData[o_idx]->addValue(SHUNT_RMIDNT, split_line[8].c_str());

    /*
     * type: real float
     * #define SHUNT_BINIT "SHUNT_BINIT"
     */
    p_busData[o_idx]->addValue(SHUNT_BINIT, atof(split_line[9].c_str()));

    if (nval > 10)
      p_busData[o_idx]->addValue(SHUNT_N1, atoi(split_line[10].c_str()));

    /*
     * type: integer
     * #define SHUNT_N2 "SHUNT_N2"
     */
    if (nval > 12)
      p_busData[o_idx]->addValue(SHUNT_N2, atoi(split_line[12].c_str()));

    /*
     * type: integer
     * #define SHUNT_N3 "SHUNT_N3"
     */
    if (nval > 14)
      p_busData[o_idx]->addValue(SHUNT_N3, atoi(split_line[14].c_str()));

    /*
     * type: integer
     * #define SHUNT_N4 "SHUNT_N4"
     */
    if (nval > 16)
      p_busData[o_idx]->addValue(SHUNT_N4, atoi(split_line[16].c_str()));

    /*
     * type: integer
     * #define SHUNT_N5 "SHUNT_N5"
     */
    if (nval > 18)
      p_busData[o_idx]->addValue(SHUNT_N5, atoi(split_line[18].c_str()));

    /*
     * type: integer
     * #define SHUNT_N6 "SHUNT_N6"
     */
    if (nval > 20)
      p_busData[o_idx]->addValue(SHUNT_N6, atoi(split_line[20].c_str()));

    /*
     * type: integer
     * #define SHUNT_N7 "SHUNT_N7"
     */
    if (nval > 22) 
      p_busData[o_idx]->addValue(SHUNT_N7, atoi(split_line[22].c_str()));

    /*
     * type: integer
     * #define SHUNT_N8 "SHUNT_N8"
     */
    if (nval > 24) 
      p_busData[o_idx]->addValue(SHUNT_N8, atoi(split_line[24].c_str()));

    /*
     * type: real float
     * #define SHUNT_B1 "SHUNT_B1"
     */
    if (nval > 11) 
      p_busData[o_idx]->addValue(SHUNT_B1, atof(split_line[11].c_str()));

    /*
     * type: real float
     * #define SHUNT_B2 "SHUNT_B2"
     */
    if (nval > 13) 
      p_busData[o_idx]->addValue(SHUNT_B2, atof(split_line[13].c_str()));

    /*
     * type: real float
     * #define SHUNT_B3 "SHUNT_B3"
     */
    if (nval > 15) 
      p_busData[o_idx]->addValue(SHUNT_B3, atof(split_line[15].c_str()));

    /*
     * type: real float
     * #define SHUNT_B4 "SHUNT_B4"
     */
    if (nval > 17) 
      p_busData[o_idx]->addValue(SHUNT_B4, atof(split_line[17].c_str()));

    /*
     * type: real float
     * #define SHUNT_B5 "SHUNT_B5"
     */
    if (nval > 19) 
      p_busData[o_idx]->addValue(SHUNT_B5, atof(split_line[19].c_str()));

    /*
     * type: real float
     * #define SHUNT_B6 "SHUNT_B6"
     */
    if (nval > 21) 
      p_busData[o_idx]->addValue(SHUNT_B6, atof(split_line[21].c_str()));

    /*
     * type: real float
     * #define SHUNT_B7 "SHUNT_B7"
     */
    if (nval > 23) 
      p_busData[o_idx]->addValue(SHUNT_B7, atof(split_line[23].c_str()));

    /*
     * type: real float
     * #define SHUNT_B8 "SHUNT_B8"
     */
    if (nval > 25) 
      p_busData[o_idx]->addValue(SHUNT_B8, atof(split_line[25].c_str()));

    stream.nextLine(line);
  }
}
