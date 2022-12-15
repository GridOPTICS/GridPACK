/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * multi_section_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "multi_section_parser33.hpp"

/**
 * Constructor
 * @param bus_map map indices in RAW file to internal indices
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to
 internal indices
 */
gridpack::parser::MultiSectParser33::MultiSectParser33(
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
gridpack::parser::MultiSectParser33::~MultiSectParser33(void)
{
}

/**
 * parse multi section block
 * @param stream input stream that feeds lines from RAW file
 * @param data vector of data collection objects to store parameters
 *             from RAW file for branches
 */
void gridpack::parser::MultiSectParser33::parse(
    gridpack::stream::InputStream &stream,
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &p_branchData)
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
    int o_idx1, o_idx2;
    o_idx1 = getBusIndex(split_line[0]);
    o_idx2 = getBusIndex(split_line[1]);

    // find branch corresponding to this line.
    int l_idx = 0;
    std::pair<int,int> branch_pair = std::pair<int,int>(o_idx1, o_idx2);
    std::map<std::pair<int, int>, int>::iterator it;
    it = p_branchMap->find(branch_pair);
    int nelems;
    bool found = false;
    if (it != p_branchMap->end()) {
      l_idx = it->second;
      found = true;
    } else {
      // Check to see if from and to buses have been switched
      std::pair<int, int> new_branch_pair;
      new_branch_pair = std::pair<int,int>(o_idx2, o_idx1);
      it = p_branchMap->find(new_branch_pair);
      if (it != p_branchMap->end()) {
        l_idx = it->second;
        found = true;
      }
    }
    nelems = split_line.size() - 3;

    if (!found || nelems <= 0) continue;

    // Clean up 2 character tag
    gridpack::utility::StringUtils util;
    std::string tag = util.clean2Char(split_line[2]);
    if (tag.length() != 2 || tag[0] != '&') {
      tag = "&1";
    }

    /*
     * type: string
     * #define MULTI_SEC_LINE_ID "MULTI_SEC_LINE_ID"
     */
    p_branchData[l_idx]->addValue(MULTI_SEC_LINE_ID, tag.c_str());

    /*
     * type: integer
     * #define MULTI_SEC_LINE_MET "MULTI_SEC_LINE_MET"
     */
    p_branchData[l_idx]->addValue(MULTI_SEC_LINE_MET, atoi(split_line[3].c_str()));


    int i;
    char buf[32];
    for (i=0; i<9; i++) {
      sprintf(buf,"MULTI_SEC_LINE_DUM%d",i+1);
      p_branchData[l_idx]->addValue(buf,atoi(split_line[i+4].c_str()));
    }

    stream.nextLine(line);
  }
}
