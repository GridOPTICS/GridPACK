// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   environment.cpp
 * @author Shrirang Abhyankar
 * @date   2022-10-05 08:52:16 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <stdio.h>
#include <ga++.h>
#include "environment.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/utilities/string_utils.hpp"
#include "gridpack/configuration/no_print.hpp"

namespace gridpack {

void Environment::PrintHelp(char** argv,const char* help)

{
  // Check if help (-h or -help) is given at command line
  if(clparser.cmdOptionExists("-h") || clparser.cmdOptionExists("-help")) {
    printf("Application Name:\n\t %s\n",argv[0]);
    if(help) printf("Description:\n\t %s\n",help);
    exit(1);
  }
}

void Environment::PrintStatus()
{
  gridpack::utility::StringUtils util;
  std::string option;
  gridpack::NoPrint* noprint = gridpack::NoPrint::instance();
  if (clparser.cmdOptionExists("-no-print")) {
    option = clparser.getCmdOption("-no-print");
  }
  if (clparser.cmdOptionExists("-noprint")) {
    option = clparser.getCmdOption("-noprint");
  }
  if (option.size() > 0) noprint->setStatus(util.getBool(option.c_str()));
}
// -------------------------------------------------------------
//  class Environment
// -------------------------------------------------------------

Environment::Environment(int argc, char **argv):p_boostEnv(argc,argv),clparser(argc,argv)
{
  pma_stack = 200000;
  pma_heap  = 200000;

  PrintStatus();
  GA_Initialize();
  MA_init(C_DBL,pma_stack,pma_heap);
  gridpack::math::Initialize(&argc,&argv);
}

Environment::Environment(int argc, char **argv,const char* help): p_boostEnv(argc,argv),clparser(argc,argv)
{
  PrintHelp(argv,help);
  pma_stack = 200000;
  pma_heap  = 200000;

  PrintStatus();
  GA_Initialize();
  MA_init(C_DBL,pma_stack,pma_heap);
  gridpack::math::Initialize(&argc,&argv);
}
  
Environment::Environment(int argc, char **argv,const char* help,const long int& ma_stack,const long int& ma_heap): p_boostEnv(argc,argv),clparser(argc,argv)
{
  PrintHelp(argv,help);

  PrintStatus();
  GA_Initialize();
  MA_init(C_DBL,ma_stack,ma_heap);
  gridpack::math::Initialize(&argc,&argv);

}

Environment::Environment(int argc, char **argv, MPI_Comm &comm): p_boostEnv(argc,argv),clparser(argc,argv)
{
  PrintStatus();
  pma_stack = 200000;
  pma_heap  = 200000;
  GA_Initialize();
  MA_init(C_DBL,pma_stack,pma_heap);
  gridpack::math::Initialize(&argc,&argv);
}

Environment::~Environment(void)
{
  // Finalize math libraries
  gridpack::math::Finalize();

  GA_Terminate();
}

  /**
   *    * return next token after option
   *       */
const std::string Environment::getCmdOption(std::string &option)
{
  return clparser.getCmdOption(option);
}

} // namespace gridpack

