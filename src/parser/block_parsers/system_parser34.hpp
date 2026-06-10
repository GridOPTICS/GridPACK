/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * system_parser34.hpp
 *       Created on: December 5, 2022
 *           Author: Bruce Palmer
 */
#ifndef _SYSTEM_PARSER34_H
#define _SYSTEM_PARSER34_H

#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/block_parsers/base_block_parser.hpp"

namespace gridpack {
namespace parser {

class SystemParser34 : public BaseBlockParser {
  public:
  SystemParser34(
      std::map<int,int> *bus_map,
      std::map<std::string,int> *name_map,
      std::map<std::pair<int, int>, int> *branch_map);

  virtual ~SystemParser34(void);

  // Parse System-Wide Data block (GENERAL, GAUSS, NEWTON, ADJUST, TYSL,
  // SOLVER, RATING, ...) and store the powerflow-relevant fields into
  // network_data. Lines that don't match a known record keyword are
  // consumed and ignored.
  void parse(
      gridpack::stream::InputStream &stream,
      boost::shared_ptr<gridpack::component::DataCollection> network_data);

  // Backward-compatible overload: just consume the block.
  void parse(
      gridpack::stream::InputStream &stream);
};

} // parser
} // gridpack
#endif
