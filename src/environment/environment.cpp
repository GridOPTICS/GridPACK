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
 * @date   2019-06-25 12:45:58 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <stdio.h>
#include "environment.hpp"
#include "gridpack/parallel/environment.hpp"
#include "gridpack/math/math.hpp"

namespace gridpack {

// -------------------------------------------------------------
//  class Environment
// -------------------------------------------------------------

// -------------------------------------------------------------
// Environment:: constructors / destructor
// -------------------------------------------------------------
  Environment::Environment(int& argc, char **argv,const char* help,const char* config_filein) : parenv(argc,argv),clparser(argc,argv)
{
  // Check if help (-h or -help) is given at command line
  if(clparser.cmdOptionExists("-h") || clparser.cmdOptionExists("-help")) {
    printf("Application Name:\n\t %s\n",argv[0]);
    if(help) printf("Description:\n\t %s\n",help);
    if(config_filein) {
      std::strcpy(config_file,config_filein);
      printf("Configuration file:\n\t%s\n",config_file);
    }
    exit(1);
  }
  gridpack::math::Initialize(&argc,&argv);
}

Environment::~Environment(void)
{
  int ierr;
  // Finalize math libraries
  gridpack::math::Finalize();

  // GA and MPI will be finalized by Environment class
}

} // namespace gridpack

