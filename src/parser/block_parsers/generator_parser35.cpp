/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * generator_parser35.cpp
 *       Created on: December 5, 2022
 *           Author: Bruce Palmer
 */
#include "generator_parser35.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to
 internal indices
 */
gridpack::parser::GeneratorParser35::GeneratorParser35(
    std::map<int,int> *bus_map,
    std::map<std::string,int> *name_map,
    std::map<std::pair<int, int>, int> *branch_map) :
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
gridpack::parser::GeneratorParser35::~GeneratorParser35(void)
{
}

/**
 * parse generator block
 * @param stream input stream that feeds lines from RAW file
 * @param data vector of data collection objects to store parameters
 *             from RAW file for generators
 */
void gridpack::parser::GeneratorParser35::parse(
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

    // GENERATOR_BUSNUMBER               "I"                   integer
    int l_idx, o_idx;
    o_idx = getBusIndex(split_line[0]);
    std::map<int, int>::iterator it;
    int nstr = split_line.size();
    it = p_busMap->find(o_idx);
    if (it != p_busMap->end()) {
      l_idx = it->second;
    } else {
      stream.nextLine(line);
      continue;
    }

    // Find out how many generators are already on bus
    int ngen;
    if (!p_busData[l_idx]->getValue(GENERATOR_NUMBER, &ngen)) ngen = 0;


    p_busData[l_idx]->addValue(GENERATOR_BUSNUMBER, o_idx, ngen);

    // Clean up 2 character tag
    gridpack::utility::StringUtils util;
    std::string tag = util.clean2Char(split_line[1]);
    // GENERATOR_ID              "ID"                  integer
    p_busData[l_idx]->addValue(GENERATOR_ID, tag.c_str(), ngen);

    // GENERATOR_PG              "PG"                  float
    if (nstr > 2) p_busData[l_idx]->addValue(GENERATOR_PG, atof(split_line[2].c_str()),
        ngen);

    // GENERATOR_QG              "QG"                  float
    if (nstr > 3) p_busData[l_idx]->addValue(GENERATOR_QG, atof(split_line[3].c_str()),
        ngen);

    // GENERATOR_QMAX              "QT"                  float
    if (nstr > 4) p_busData[l_idx]->addValue(GENERATOR_QMAX,
        atof(split_line[4].c_str()), ngen);

    // GENERATOR_QMIN              "QB"                  float
    if (nstr > 5) p_busData[l_idx]->addValue(GENERATOR_QMIN,
        atof(split_line[5].c_str()), ngen);

    // GENERATOR_VS              "VS"                  float
    if (nstr > 6) p_busData[l_idx]->addValue(GENERATOR_VS, atof(split_line[6].c_str()),
        ngen);

    // GENERATOR_IREG            "IREG"                integer
    if (nstr > 7) p_busData[l_idx]->addValue(GENERATOR_IREG,
        atoi(split_line[7].c_str()), ngen);

    // GENERATOR_NREG            "NREG"                integer
    if (nstr > 8) p_busData[l_idx]->addValue(GENERATOR_NREG,
        atoi(split_line[8].c_str()), ngen);

    // GENERATOR_MBASE           "MBASE"               float
    if (nstr > 9) p_busData[l_idx]->addValue(GENERATOR_MBASE,
        atof(split_line[9].c_str()), ngen);

    // GENERATOR_ZSOURCE                                complex
    if (nstr > 11) p_busData[l_idx]->addValue(GENERATOR_ZSOURCE,
        gridpack::ComplexType(atof(split_line[10].c_str()),
          atof(split_line[11].c_str())), ngen);

    // GENERATOR_XTRAN                              complex
    if (nstr > 13) p_busData[l_idx]->addValue(GENERATOR_XTRAN,
        gridpack::ComplexType(atof(split_line[12].c_str()),
          atof(split_line[13].c_str())), ngen);

    // GENERATOR_GTAP              "GTAP"                  float
    if (nstr > 14) p_busData[l_idx]->addValue(GENERATOR_GTAP,
        atof(split_line[14].c_str()), ngen);

    // GENERATOR_STAT              "STAT"                  float
    if (nstr > 15)  p_busData[l_idx]->addValue(GENERATOR_STAT,
        atoi(split_line[15].c_str()), ngen);

    // GENERATOR_RMPCT           "RMPCT"               float
    if (nstr > 16) p_busData[l_idx]->addValue(GENERATOR_RMPCT,
        atof(split_line[16].c_str()), ngen);

    // GENERATOR_PMAX              "PT"                  float
    if (nstr > 17) p_busData[l_idx]->addValue(GENERATOR_PMAX,
        atof(split_line[17].c_str()), ngen);

    // GENERATOR_PMIN              "PB"                  float
    if (nstr > 18) p_busData[l_idx]->addValue(GENERATOR_PMIN,
        atof(split_line[18].c_str()), ngen);

    // GENERATOR_BASLOD              "BASLOD"            integer 
    if (nstr > 19) p_busData[l_idx]->addValue(GENERATOR_BASLOD,
        atof(split_line[19].c_str()), ngen);

    // TODO: add variables Oi, Fi, WMOD, WPF
    // There may be between 0 and 4 owner pairs.
    // In order for WMOD and WPF to be defined, there need to be 28 entries
    // in the line. The owners may be included as blanks.
    if (nstr > 20) {
      int owner = 0;
      if (this->isBlank(split_line[20])) {
        p_busData[l_idx]->getValue(BUS_OWNER,&owner);
      } else {
        owner = atoi(split_line[20].c_str());
      }
      p_busData[l_idx]->addValue(GENERATOR_OWNER1, owner, ngen);
      double frac = 1.0;
      if (nstr > 21) {
        if (!this->isBlank(split_line[21])) {
          frac = atof(split_line[21].c_str());
        }
      }
      p_busData[l_idx]->addValue(GENERATOR_OFRAC1, frac, ngen);
    }
    if (nstr > 22) {
      int owner = 0;
      if (this->isBlank(split_line[22])) {
        owner = 0;
      } else {
        owner = atoi(split_line[22].c_str());
      }
      p_busData[l_idx]->addValue(GENERATOR_OWNER2, owner, ngen);
      double frac = 0.0;
      if (nstr > 23) {
        if (!this->isBlank(split_line[23])) {
          frac = atof(split_line[23].c_str());
        }
      }
      p_busData[l_idx]->addValue(GENERATOR_OFRAC2, frac, ngen);
    }
    if (nstr > 24) {
      int owner = 0;
      if (this->isBlank(split_line[24])) {
        owner = 0;
      } else {
        owner = atoi(split_line[24].c_str());
      }
      p_busData[l_idx]->addValue(GENERATOR_OWNER3, owner, ngen);
      double frac = 0.0;
      if (nstr > 25) {
        if (!this->isBlank(split_line[25])) {
          frac = atof(split_line[25].c_str());
        }
      }
      p_busData[l_idx]->addValue(GENERATOR_OFRAC3, frac, ngen);
    }
    if (nstr > 26) {
      int owner = 0;
      if (this->isBlank(split_line[26])) {
        owner = 0;
      } else {
        owner = atoi(split_line[26].c_str());
      }
      p_busData[l_idx]->addValue(GENERATOR_OWNER4, owner, ngen);
      double frac = 0.0;
      if (nstr > 27) {
        if (!this->isBlank(split_line[27])) {
          frac = atof(split_line[27].c_str());
        }
      }
      p_busData[l_idx]->addValue(GENERATOR_OFRAC4, frac, ngen);
    }

    // Last two entries are WMOD and WPF
    if (nstr > 28) {
      p_busData[l_idx]->addValue(GENERATOR_WMOD,
          atoi(split_line[28].c_str()), ngen);
    }
    if (nstr > 29) {
      p_busData[l_idx]->addValue(GENERATOR_WPF,
          atof(split_line[29].c_str()), ngen);
    }

    // Increment number of generators in data object
    if (ngen == 0) {
      ngen = 1;
      p_busData[l_idx]->addValue(GENERATOR_NUMBER,ngen);
    } else {
      ngen++;
      p_busData[l_idx]->setValue(GENERATOR_NUMBER,ngen);
    }

    stream.nextLine(line);
  }
}
