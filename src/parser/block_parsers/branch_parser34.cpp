/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * branch_parser34.cpp
 *       Created on: December 5, 2022
 *           Author: Bruce Palmer
 */
#include "branch_parser34.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to
 internal indices
 */
gridpack::parser::BranchParser34::BranchParser34(
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
gridpack::parser::BranchParser34::~BranchParser34(void)
{
}

/**
 * parse branch block
 * @param stream input stream that feeds lines from RAW file
 * @param data vector of data collection objects to store parameters
 *             from RAW file for branches
 */
void gridpack::parser::BranchParser34::parse(
    gridpack::stream::InputStream &stream,
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &p_branchData)
{
  std::string line;
  int  o_idx1, o_idx2;
  int index = 0;

  stream.nextLine(line); //this should be the first line of the block

  int nelems;
  while(test_end(line)) {
    std::pair<int, int> branch_pair;
    std::vector<std::string>  split_line;
    if (check_comment(line)) {
      stream.nextLine(line);
      continue;
    }
    this->cleanComment(line);
    split_line = this->splitPSSELine(line);

    o_idx1 = getBusIndex(split_line[0]);
    o_idx2 = getBusIndex(split_line[1]);

    // Check to see if pair has already been created
    int l_idx = 0;
    branch_pair = std::pair<int,int>(o_idx1, o_idx2);
    std::map<std::pair<int, int>, int>::iterator it;
    it = p_branchMap->find(branch_pair);

    bool switched = false;
    if (it != p_branchMap->end()) {
      l_idx = it->second;
      p_branchData[l_idx]->getValue(BRANCH_NUM_ELEMENTS,&nelems);
    } else {
      // Check to see if from and to buses have been switched
      std::pair<int, int> new_branch_pair;
      new_branch_pair = std::pair<int,int>(o_idx2, o_idx1);
      it = p_branchMap->find(new_branch_pair);
      if (it != p_branchMap->end()) {
        //            printf("Found multiple lines with switched buses 1: %d 2: %d\n",
        //                o_idx1,o_idx2);
        l_idx = it->second;
        p_branchData[l_idx]->getValue(BRANCH_NUM_ELEMENTS,&nelems);
        switched = true;
      } else {
        boost::shared_ptr<gridpack::component::DataCollection>
          data(new gridpack::component::DataCollection);
        l_idx = p_branchData.size();
        p_branchData.push_back(data);
        nelems = 0;
        p_branchData[l_idx]->addValue(BRANCH_NUM_ELEMENTS,nelems);
      }
    }

    if (nelems == 0) {
      // BRANCH_INDEX                                   integer
      p_branchData[l_idx]->addValue(BRANCH_INDEX, index);

      // BRANCH_FROMBUS            "I"                   integer
      p_branchData[l_idx]->addValue(BRANCH_FROMBUS, o_idx1);

      // BRANCH_TOBUS            "J"                   integer
      p_branchData[l_idx]->addValue(BRANCH_TOBUS, o_idx2);

      // add pair to branch map
      p_branchMap->insert(std::pair<std::pair<int, int>, int >(branch_pair,
            index));
      index++;
    }

    // BRANCH_SWITCHED
    p_branchData[l_idx]->addValue(BRANCH_SWITCHED, switched, nelems);

    int nstr = split_line.size();
    // Clean up 2 character tag
    gridpack::utility::StringUtils util;
    std::string tag = util.clean2Char(split_line[2]);
    // BRANCH_CKT          "CKT"                 character
    if (nstr > 2) p_branchData[l_idx]->addValue(BRANCH_CKT,
        tag.c_str(), nelems);

    // BRANCH_R            "R"                   float
    if (nstr > 3) p_branchData[l_idx]->addValue(BRANCH_R,
        atof(split_line[3].c_str()), nelems);

    // BRANCH_X            "X"                   float
    if (nstr > 4) p_branchData[l_idx]->addValue(BRANCH_X,
        atof(split_line[4].c_str()), nelems);

    // BRANCH_B            "B"                   float
    if (nstr > 5) p_branchData[l_idx]->addValue(BRANCH_B,
        atof(split_line[5].c_str()), nelems);


    // BRANCH_NAME             "NAME"                 string
    if (nstr > 6) {
      if (split_line[6].find_first_of('\'',0) != std::string::npos) {
        gridpack::utility::StringUtils util;
        util.trim(split_line[6]);
        p_branchData[l_idx]->addValue(BRANCH_NAME, split_line[6].c_str(), nelems);
      }
    }
    // BRANCH_RATE1-12                              float
    if (nstr > 7) {
      p_branchData[l_idx]->addValue(BRANCH_RATE1,
        atof(split_line[7].c_str()), nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATING_A,
        atof(split_line[7].c_str()), nelems);
    }
    if (nstr > 8) {
      p_branchData[l_idx]->addValue(BRANCH_RATE2,
        atof(split_line[8].c_str()), nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATING_B,
        atof(split_line[8].c_str()), nelems);
    }
    if (nstr > 9) {
      p_branchData[l_idx]->addValue(BRANCH_RATE3,
        atof(split_line[9].c_str()), nelems);
      p_branchData[l_idx]->addValue(BRANCH_RATING_C,
        atof(split_line[9].c_str()), nelems);
    }
    if (nstr > 10) p_branchData[l_idx]->addValue(BRANCH_RATE4,
        atof(split_line[10].c_str()), nelems);
    if (nstr > 11) p_branchData[l_idx]->addValue(BRANCH_RATE5,
        atof(split_line[11].c_str()), nelems);
    if (nstr > 12) p_branchData[l_idx]->addValue(BRANCH_RATE6,
        atof(split_line[12].c_str()), nelems);
    if (nstr > 13) p_branchData[l_idx]->addValue(BRANCH_RATE7,
        atof(split_line[13].c_str()), nelems);
    if (nstr > 14) p_branchData[l_idx]->addValue(BRANCH_RATE8,
        atof(split_line[14].c_str()), nelems);
    if (nstr > 15) p_branchData[l_idx]->addValue(BRANCH_RATE9,
        atof(split_line[15].c_str()), nelems);
    if (nstr > 16) p_branchData[l_idx]->addValue(BRANCH_RATE10,
        atof(split_line[16].c_str()), nelems);
    if (nstr > 17) p_branchData[l_idx]->addValue(BRANCH_RATE11,
        atof(split_line[17].c_str()), nelems);
    if (nstr > 18) p_branchData[l_idx]->addValue(BRANCH_RATE12,
        atof(split_line[18].c_str()), nelems);

    // BRANCH_SHUNT_ADMTTNC_G1        "GI"               float
    if (nstr > 19) p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_G1,
        atof(split_line[19].c_str()), nelems);

    // BRANCH_SHUNT_ADMTTNC_B1        "BI"               float
    if (nstr > 20) p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_B1,
        atof(split_line[20].c_str()), nelems);

    // BRANCH_SHUNT_ADMTTNC_G2        "GJ"               float
    if (nstr > 21) p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_G2,
        atof(split_line[21].c_str()), nelems);

    // BRANCH_SHUNT_ADMTTNC_B2        "BJ"               float
    if (nstr > 22) p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_B2,
        atof(split_line[22].c_str()), nelems);

    // BRANCH_STATUS        "STATUS"               integer
    if (nstr > 23) p_branchData[l_idx]->addValue(BRANCH_STATUS,
        atoi(split_line[23].c_str()), nelems);

    // BRANCH_METER         "MET"                  integer
    if (nstr > 24) p_branchData[l_idx]->addValue(BRANCH_METER,
        atoi(split_line[24].c_str()), nelems);

    // BRANCH_LENGTH        "LEN"                        float
    if (nstr > 25) p_branchData[l_idx]->addValue(BRANCH_LENGTH,
        atof(split_line[25].c_str()), nelems);

    // BRANCH_O1        "O1"                       integer
    if (nstr > 26) p_branchData[l_idx]->addValue(BRANCH_O1,
        atoi(split_line[26].c_str()), nelems);

    // BRANCH_F1        "F1"                             float
    if (nstr > 27) p_branchData[l_idx]->addValue(BRANCH_F1,
        atoi(split_line[27].c_str()), nelems);

    // BRANCH_O2        "O2"                       integer
    if (nstr > 28) p_branchData[l_idx]->addValue(BRANCH_O2,
        atoi(split_line[28].c_str()), nelems);

    // BRANCH_F2        "F2"                             float
    if (nstr > 29) p_branchData[l_idx]->addValue(BRANCH_F2,
        atoi(split_line[29].c_str()), nelems);

    // BRANCH_O3        "O3"                       integer
    if (nstr > 30) p_branchData[l_idx]->addValue(BRANCH_O3,
        atoi(split_line[30].c_str()), nelems);

    // BRANCH_F3        "F3"                             float
    if (nstr > 31) p_branchData[l_idx]->addValue(BRANCH_F3,
        atoi(split_line[31].c_str()), nelems);

    // BRANCH_O4        "O4"                       integer
    if (nstr > 32) p_branchData[l_idx]->addValue(BRANCH_O4,
        atoi(split_line[32].c_str()), nelems);

    // BRANCH_F4        "F4"                             float
    if (nstr > 33) p_branchData[l_idx]->addValue(BRANCH_F4,
        atoi(split_line[33].c_str()), nelems);

    // TODO: add variables MET, LEN, Oi, Fi

    // Add BRANCH_TAP with value 0.0
    p_branchData[l_idx]->addValue(BRANCH_TAP,0.0,nelems);

    // Add BRANCH_SHIFT with value 0.0
    p_branchData[l_idx]->addValue(BRANCH_SHIFT,0.0,nelems);

    nelems++;
    p_branchData[l_idx]->setValue(BRANCH_NUM_ELEMENTS,nelems);
    stream.nextLine(line);
  }
}
