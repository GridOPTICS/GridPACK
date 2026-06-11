/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 *
 *
 * two_term_parser33.hpp
 *       Created on: November 29, 2022
 *           Author: Bruce Palmer
 */
#ifndef _TWO_TERM_PARSER33_H
#define _TWO_TERM_PARSER33_H

#include "gridpack/parser/block_parsers/base_block_parser.hpp"

namespace gridpack {
namespace parser {

class TwoTermParser33 : public BaseBlockParser {
  public:
  /**
   * Constructor
   * @param bus_map map indices in RAW file to internal indices
   * @param name_map map name in RAW file to internal indices
   * @param branch_map map bus index pair in RAW file to internal indices
   */
  TwoTermParser33(
      std::map<int,int> *bus_map,
      std::map<std::string,int> *name_map,
      std::map<std::pair<int, int>, int> *branch_map);

  /**
   * Simple Destructor
   */
  virtual ~TwoTermParser33(void);

  /**
   * parse two terminal block. Currently does not store data.
   * Kept for callers that don't need the DC injection synthesized into loads.
   * @param stream input stream that feeds lines from RAW file
   */
  void parse(
      gridpack::stream::InputStream &stream);

  /**
   * Parse two-terminal DC block and synthesize PQ injections at rectifier and
   * inverter as additional LOAD records on the corresponding buses. PSS/E
   * power-flow init holds DC at scheduled MW; if we discard the records the
   * AC system sees a phantom mismatch of |SETVL| at each terminal and nearby
   * generators are driven to Q limits.
   *
   * Modeling:
   *  - Rectifier bus: positive P load = SETVL (rectifier consumes from AC)
   *    plus reactive load = 0.5 * SETVL (rule-of-thumb 50% of P for converter VARs)
   *  - Inverter bus: negative P load = -SETVL (inverter supplies to AC)
   *    plus reactive load = 0.5 * SETVL
   *  - Records with MDC = 0 (blocked) are skipped.
   *  - SETVL is in MW (when MDC=1, current control) or kA (when MDC=2, current).
   *    For MDC=2 we estimate P = SETVL_kA * VSCHD_kV as MW.
   *
   * @param stream input stream
   * @param busData per-bus DataCollection vector to which LOAD records are appended
   * @param case_sbase ignored (loads are stored in MW/MVAr per LOAD format)
   */
  void parse(
      gridpack::stream::InputStream &stream,
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> > &busData,
      double case_sbase);
};

} // parser
} // gridpack
#endif
