/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * case_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "case_parser33.hpp"

/**
 * Constructor
 * @param stream input stream that feeds lines from RAW file
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to
 internal indices
 */
gridpack::parser::CaseParser33::CaseParser33(
    std::map<int,int> &bus_map,
    std::map<std::string,int> &name_map,
    std::map<std::pair<int, int>, int> &branch_map) :
    gridpack::parser::BaseBlockParser(
      bus_map, name_map, branch_map)
{
}


/**
 * Simple Destructor
 */
gridpack::parser::CaseParser33::~CaseParser33(void)
{
}

/**
 * parse case block
 * @param bus_map map indices in RAW file to internal indices
 * @param data data collection object to store parameters from RAW file
 * @param sbase value of sbase from RAW file
 * @param id value if id from RAW file
 */
void gridpack::parser::CaseParser33::parse(
    gridpack::stream::InputStream &stream,
    boost::shared_ptr<gridpack::component::DataCollection> &data,
    double &sbase, int &id)
{
  std::string line;

  stream.nextLine(line);
  while (check_comment(line)) {
    stream.nextLine(line);
  }
  std::vector<std::string>  split_line;

  this->cleanComment(line);
  boost::algorithm::split(split_line, line, boost::algorithm::is_any_of(","),
      boost::token_compress_off);

  // CASE_ID             "IC"                   ranged integer
  id = atoi(split_line[0].c_str());

  // CASE_SBASE          "SBASE"                float
  sbase = atof(split_line[1].c_str());

  data->addValue(CASE_SBASE, p_case_sbase);
  data->addValue(CASE_ID, p_case_id);
  /*  These do not appear in the dictionary
  // REVISION_ID
  if (split_line.size() > 2)
  p_revision_id = atoi(split_line[2].c_str());

  // XFRRAT_UNITS
  if (split_line.size() > 3)
  p_xffrat_units = atof(split_line[3].c_str());

  // NXFRAT_UNITS
  if (split_line.size() > 4)
  p_nxfrat_units = atof(split_line[4].c_str());

  // BASE_FREQ
  if (split_line.size() > 5)
  p_base_freq = atof(split_line[5].c_str());
  */

}
