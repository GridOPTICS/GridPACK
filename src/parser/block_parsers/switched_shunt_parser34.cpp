/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * switched_shunt_parser34.cpp
 *       Created on: December 9, 2022
 *           Author: Bruce Palmer
 */
#include "switched_shunt_parser34.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to
 internal indices
 */
gridpack::parser::SwitchedShuntParser34::SwitchedShuntParser34(
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
gridpack::parser::SwitchedShuntParser34::~SwitchedShuntParser34(void)
{
}

/**
 * parse switched shunt block
 * @param stream input stream that feeds lines from RAW file
 * @param data vector of data collection objects to store parameters
 *             from RAW file for switched shunts
 */
void gridpack::parser::SwitchedShuntParser34::parse(
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
    std::map<int, int>::iterator it;
    it = p_busMap->find(l_idx);
    if (it != p_busMap->end()) {
      o_idx = it->second;
    } else {
      stream.nextLine(line);
      continue;
    }
    int nval = split_line.size();

    if (!p_busData[o_idx]->setValue(SWSHUNT_BUSNUMBER,
          atoi(split_line[0].c_str()))) {
      p_busData[o_idx]->addValue(SWSHUNT_BUSNUMBER, atoi(split_line[0].c_str()));
    }

    // Currently ignoring the shunt ID. Assume one switched shunt per bus
    /*
     * type: integer
     * #define SHUNT_MODSW "SHUNT_MODSW"
     */
    if (!p_busData[o_idx]->setValue(SHUNT_MODSW,
          atoi(split_line[1].c_str()))) {
      p_busData[o_idx]->addValue(SHUNT_MODSW, atoi(split_line[1].c_str()));
    }

    /*
     * type: integer
     * #define SHUNT_ADJM "SHUNT_ADJM"
     */
    if (!p_busData[o_idx]->setValue(SHUNT_ADJM,
          atoi(split_line[2].c_str()))) {
      p_busData[o_idx]->addValue(SHUNT_ADJM, atoi(split_line[2].c_str()));
    }

    /*
     * type: integer
     * #define SHUNT_SWCH_STAT "SHUNT_SWCH_STAT"
     */
    if (!p_busData[o_idx]->setValue(SHUNT_SWCH_STAT,
          atoi(split_line[3].c_str()))) {
      p_busData[o_idx]->addValue(SHUNT_SWCH_STAT, atoi(split_line[3].c_str()));
    }

    /*
     * type: real float
     * #define SHUNT_VSWHI "SHUNT_VSWHI"
     */
    if (!p_busData[o_idx]->setValue(SHUNT_VSWHI,
          atof(split_line[4].c_str()))) {
      p_busData[o_idx]->addValue(SHUNT_VSWHI, atof(split_line[4].c_str()));
    }

    /*
     * type: real float
     * #define SHUNT_VSWLO "SHUNT_VSWLO"
     */
    if (!p_busData[o_idx]->setValue(SHUNT_VSWLO,
          atof(split_line[5].c_str()))) {
      p_busData[o_idx]->addValue(SHUNT_VSWLO, atof(split_line[5].c_str()));
    }

    /*
     * type: integer
     * #define SHUNT_SWREG "SHUNT_SWREG"
     */
    if (!p_busData[o_idx]->setValue(SHUNT_SWREG,
          atoi(split_line[6].c_str()))) {
      p_busData[o_idx]->addValue(SHUNT_SWREG, atoi(split_line[6].c_str()));
    }

    /*
     * type: real float
     * #define SHUNT_RMPCT "SHUNT_RMPCT"
     */
    if (!p_busData[o_idx]->setValue(SHUNT_RMPCT,
          atof(split_line[7].c_str()))) {
      p_busData[o_idx]->addValue(SHUNT_RMPCT, atof(split_line[7].c_str()));
    }

    /*
     * type: string
     * #define SHUNT_RMIDNT "SHUNT_RMIDNT"
     */
    if (!p_busData[o_idx]->setValue(SHUNT_RMIDNT,
          split_line[8].c_str())) {
      p_busData[o_idx]->addValue(SHUNT_RMIDNT, split_line[8].c_str());
    }

    /*
     * type: real float
     * #define SHUNT_BINIT "SHUNT_BINIT"
     */
    if (!p_busData[o_idx]->setValue(SHUNT_BINIT,
          atof(split_line[9].c_str()))) {
      p_busData[o_idx]->addValue(SHUNT_BINIT, atof(split_line[9].c_str()));
    }

    /*
     * type: integer
     * #define SHUNT_N1 "SHUNT_N1"
     */
    if (nval > 10) {
      if (!p_busData[o_idx]->setValue(SHUNT_N1,
            atoi(split_line[10].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_N1, atoi(split_line[10].c_str()));
      }
    }

    /*
     * type: integer
     * #define SHUNT_N2 "SHUNT_N2"
     */
    if (nval > 12) {
      if (!p_busData[o_idx]->setValue(SHUNT_N2,
            atoi(split_line[12].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_N2, atoi(split_line[12].c_str()));
      }
    }

    /*
     * type: integer
     * #define SHUNT_N3 "SHUNT_N3"
     */
    if (nval > 14) {
      if (!p_busData[o_idx]->setValue(SHUNT_N3,
            atoi(split_line[14].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_N3, atoi(split_line[14].c_str()));
      }
    }

    /*
     * type: integer
     * #define SHUNT_N4 "SHUNT_N4"
     */
    if (nval > 16) {
      if (!p_busData[o_idx]->setValue(SHUNT_N4,
            atoi(split_line[16].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_N4, atoi(split_line[16].c_str()));
      }
    }

    /*
     * type: integer
     * #define SHUNT_N5 "SHUNT_N5"
     */
    if (nval > 18) {
      if (!p_busData[o_idx]->setValue(SHUNT_N5,
            atoi(split_line[18].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_N5, atoi(split_line[18].c_str()));
      }
    }

    /*
     * type: integer
     * #define SHUNT_N6 "SHUNT_N6"
     */
    if (nval > 20) {
      if (!p_busData[o_idx]->setValue(SHUNT_N6,
            atoi(split_line[20].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_N6, atoi(split_line[20].c_str()));
      }
    }

    /*
     * type: integer
     * #define SHUNT_N7 "SHUNT_N7"
     */
    if (nval > 22) {
      if (!p_busData[o_idx]->setValue(SHUNT_N7,
            atoi(split_line[22].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_N7, atoi(split_line[22].c_str()));
      }
    }

    /*
     * type: integer
     * #define SHUNT_N8 "SHUNT_N8"
     */
    if (nval > 24) {
      if (!p_busData[o_idx]->setValue(SHUNT_N8,
            atoi(split_line[24].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_N8, atoi(split_line[24].c_str()));
      }
    }

    /*
     * type: real float
     * #define SHUNT_B1 "SHUNT_B1"
     */
    if (nval > 11) {
      if (!p_busData[o_idx]->setValue(SHUNT_B1,
            atof(split_line[11].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_B1, atof(split_line[11].c_str()));
      }
    }

    /*
     * type: real float
     * #define SHUNT_B2 "SHUNT_B2"
     */
    if (nval > 13) {
      if (!p_busData[o_idx]->setValue(SHUNT_B2,
            atof(split_line[13].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_B2, atof(split_line[13].c_str()));
      }
    }

    /*
     * type: real float
     * #define SHUNT_B3 "SHUNT_B3"
     */
    if (nval > 15) {
      if (!p_busData[o_idx]->setValue(SHUNT_B3,
            atof(split_line[15].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_B3, atof(split_line[15].c_str()));
      }
    }

    /*
     * type: real float
     * #define SHUNT_B4 "SHUNT_B4"
     */
    if (nval > 17) {
      if (!p_busData[o_idx]->setValue(SHUNT_B4,
            atof(split_line[17].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_B4, atof(split_line[17].c_str()));
      }
    }


    /*
     * type: real float
     * #define SHUNT_B5 "SHUNT_B5"
     */
    if (nval > 19) {
      if (!p_busData[o_idx]->setValue(SHUNT_B5,
            atof(split_line[19].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_B5, atof(split_line[19].c_str()));
      }
    }

    /*
     * type: real float
     * #define SHUNT_B6 "SHUNT_B6"
     */
    if (nval > 21) {
      if (!p_busData[o_idx]->setValue(SHUNT_B6,
            atof(split_line[21].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_B6, atof(split_line[21].c_str()));
      }
    }

    /*
     * type: real float
     * #define SHUNT_B7 "SHUNT_B7"
     */
    if (nval > 23) {
      if (!p_busData[o_idx]->setValue(SHUNT_B7,
            atof(split_line[23].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_B7, atof(split_line[23].c_str()));
      }
    }

    /*
     * type: real float
     * #define SHUNT_B8 "SHUNT_B8"
     */
    if (nval > 25) {
      if (!p_busData[o_idx]->setValue(SHUNT_B8,
            atof(split_line[25].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_B8, atof(split_line[25].c_str()));
      }
    }

    /*
     * type: integer
     * #define SHUNT_NREG"
     */
    if (nval > 26) {
      if (!p_busData[o_idx]->setValue(SHUNT_NREG,
            atoi(split_line[26].c_str()))) {
        p_busData[o_idx]->addValue(SHUNT_NREG, atof(split_line[26].c_str()));
      }
    }

    stream.nextLine(line);
  }
}
