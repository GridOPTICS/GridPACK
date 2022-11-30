/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * bus_parser33.cpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#include "bus_parser33.hpp"

/**
 * Constructor
 * @param stream input stream that feeds lines from RAW file
 * @param name_map map name in RAW file to internal indices
 * @param branch_map map bus index pair in RAW file to
 internal indices
 */
gridpack::parser::BusParser33::BusParser33(
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
gridpack::parser::BusParser33::~BusParser33(void)
{
}

/**
 * parse bus block
 * @param bus_map map indices in RAW file to internal indices
 * @param data vector of data collection objects to store parameters
 *             from RAW file for buses
 * @param p_case_sbase value of sbase from RAW file
 * @param p_case_id value if id from RAW file
 * @param p_maxBusIndex maximum value of the bus indices
 */
void gridpack::parser::BusParser33::parse(
    gridpack::stream::InputStream &stream,
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &p_busData,
    double p_case_sbase, int p_case_id, int *p_maxBusIndex)
{
  std::string          line;
  int                  index = 0;
  int                  o_idx;
  stream.nextLine(line);
  stream.nextLine(line);
  stream.nextLine(line);

  while(test_end(line)) {
    std::vector<std::string>  split_line;
    if (check_comment(line)) {
      stream.nextLine(line);
      continue;
    }
    this->cleanComment(line);
    split_line = this->splitPSSELine(line);
    boost::shared_ptr<gridpack::component::DataCollection>
      data(new gridpack::component::DataCollection);
    int nstr = split_line.size();

    // BUS_I               "I"                   integer
    o_idx = atoi(split_line[0].c_str());
    if (*p_maxBusIndex<o_idx) *p_maxBusIndex = o_idx;
    data->addValue(BUS_NUMBER, o_idx);
    p_busData.push_back(data);
    p_busMap.insert(std::pair<int,int>(o_idx,index));

    // Add case parameters to data set
    data->addValue(CASE_SBASE, p_case_sbase);
    data->addValue(CASE_ID, p_case_id);

    // BUS_NAME             "NAME"                 string
    std::string bus_name = split_line[1];

    //store bus and index as a pair
    {
      check_string(bus_name);
      std::pair<std::string,int> name_pair;
      name_pair = std::pair<std::string,int>(bus_name,abs(o_idx));
      p_nameMap.insert(name_pair);
    }
    if (split_line[1].find_first_of('\'',0) != std::string::npos) {
      gridpack::utility::StringUtils util;
      util.trim(split_line[1]);
    }
    if (nstr > 1) data->addValue(BUS_NAME, bus_name.c_str());

    // BUS_BASEKV           "BASKV"               float
    if (nstr > 2) data->addValue(BUS_BASEKV, atof(split_line[2].c_str()));

    // BUS_TYPE               "IDE"                   integer
    if (nstr > 3) data->addValue(BUS_TYPE, atoi(split_line[3].c_str()));

    // BUS_AREA            "IA"                integer
    if (nstr > 4) data->addValue(BUS_AREA, atoi(split_line[4].c_str()));

    // BUS_ZONE            "ZONE"                integer
    if (nstr > 5) data->addValue(BUS_ZONE, atoi(split_line[5].c_str()));

    // BUS_OWNER              "IA"                  integer
    if (nstr > 6) data->addValue(BUS_OWNER, atoi(split_line[6].c_str()));

    // BUS_VOLTAGE_MAG              "VM"                  float
    if (nstr > 7) data->addValue(BUS_VOLTAGE_MAG, atof(split_line[7].c_str()));

    // BUS_VOLTAGE_ANG              "VA"                  float
    if (nstr > 8) data->addValue(BUS_VOLTAGE_ANG, atof(split_line[8].c_str()));

    // BUS_VOLTAGE_MAX              "VOLTAGE_MAX"               float
    if (nstr > 9) data->addValue(BUS_VOLTAGE_MAX, atof(split_line[9].c_str()));

    // BUS_VOLTAGE_MIN              "VOLTAGE_MIN"              float
    if (nstr > 10) data->addValue(BUS_VOLTAGE_MIN, atof(split_line[10].c_str()));

    // TODO: Need to add EVHI, EVLO
    index++;
    stream.nextLine(line);
  }
}
