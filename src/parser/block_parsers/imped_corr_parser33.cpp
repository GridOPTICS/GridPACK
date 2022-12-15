/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * imped_corr_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "imped_corr_parser33.hpp"

/**
 * Constructor
 * @param stream input stream that feeds lines from RAW file
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to
 internal indices
 */
gridpack::parser::ImpedCorrParser33::ImpedCorrParser33(
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
gridpack::parser::ImpedCorrParser33::~ImpedCorrParser33(void)
{
}

/**
 * parse impedence correction block
 * @param bus_map map indices in RAW file to internal indices
 * @param data vector of data collection objects to store parameters
 *             from RAW file for impedence correction table
 */
void gridpack::parser::ImpedCorrParser33::parse(
    gridpack::stream::InputStream &stream,
    std::map<int,boost::shared_ptr<gridpack::component::DataCollection> > &p_imp_corr_table)
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
    int nval = split_line.size();
    int entries = nval-1;
    entries =  entries - entries%2;
    entries = entries/2;
    if (entries == 0) continue;
    boost::shared_ptr<gridpack::component::DataCollection>
      data(new gridpack::component::DataCollection);

    /*
     * type: integer
     * #define XFMR_CORR_TABLE_NUMBER "XFMR_CORR_TABLE_NUMBER"
     */
    int tableid = atoi(split_line[0].c_str());
    data->addValue(XFMR_CORR_TABLE_NUMBER, tableid);

    int i;
    char buf[32];
    for (i=0; i<entries; i++) {
      /*
       * type: real float
       * #define XFMR_CORR_TABLE_Ti "XFMR_CORR_TABLE_Ti"
       */
      sprintf(buf,"XFMR_CORR_TABLE_T%d",i+1);
      data->addValue(buf, atof(split_line[1+2*i].c_str()));

      /*
       * type: real float
       * #define XFMR_CORR_TABLE_Fi "XFMR_CORR_TABLE_Fi"
       */
      sprintf(buf,"XFMR_CORR_TABLE_F%d",i+1);
      data->addValue(buf, atof(split_line[2+2*i].c_str()));
    }

    p_imp_corr_table.insert(std::pair<int,
        boost::shared_ptr<gridpack::component::DataCollection> >(tableid,data));
    stream.nextLine(line);
  }
}
