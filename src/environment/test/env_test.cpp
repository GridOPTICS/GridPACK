
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   env_test.cpp
 * @author Bruce Palmer
 * @date   2016-07-14 14:23:07 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "gridpack/include/gridpack.hpp"

const char* help = "Test application for Environment class. To test command\n"
                   "         line options, run this test with options -help and then\n"
                   "         run it with options -i input -o output";

int main(int argc, char **argv)
{
  // Initialize libraries (parallel and math)
  gridpack::Environment env(argc,argv,help);

  std::string option = "-i";
  std::string i_option = env.getCmdOption(option);
  std::cout<< "Option value for -i is: "<<i_option<<std::endl;
  option = "-o";
  std::string o_option = env.getCmdOption(option);
  std::cout<< "Option value for -o is: "<<o_option<<std::endl;

  return 0;
}

