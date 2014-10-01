/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   parser_c.cpp
 * @author Bruce Palmer
 * @date   2014-09-17 11:05:08 d3g293
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/parser/PTI23_parser.hpp"
#include "../component/fortran_component.hpp"

typedef gridpack::network::BaseNetwork<
  gridpack::fortran_component::FortranBusComponent,
  gridpack::fortran_component::FortranBranchComponent>
  FortranNetwork;

typedef gridpack::parser::PTI23_parser<FortranNetwork> FortranPTI23Parser;

struct networkWrapper {
    boost::shared_ptr<FortranNetwork> network;
};

/**
 * Create a PTI23 parser
 * @param parser pointer to Fortran PTI23 parser object
 * @param network pointer to Fortran network object
 */
extern "C" void pti23_parser_create(FortranPTI23Parser **parser,
    networkWrapper *wnetwork)
{
  *parser = new FortranPTI23Parser(wnetwork->network);
}

/**
 * Destroy a PTI23 parser
 * @param parser pointer to Fortran PTI23 parser object
 */
extern "C" void pti23_parser_destroy(FortranPTI23Parser **parser)
{
  delete (*parser);
}

/**
 * Parse network configuration file and create a network
 * @param parser pointer to Fortran PTI23 parser object
 * @param string name of network file
 */
extern "C" void pti23_parser_parse(FortranPTI23Parser *parser,
    char *string)
{
  parser->parse(string);
}
