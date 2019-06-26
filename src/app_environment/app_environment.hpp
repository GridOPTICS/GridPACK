// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   app_env.hpp
 * @author Shrirang Abhyankar
 * @date   2019-06-25 14:06:28 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _app_environment_hpp_
#define _app_environment_hpp_

#include "gridpack/parallel/environment.hpp"

namespace gridpack {
namespace app_environment {

class CommandLineParser {
public:
  CommandLineParser (int &argc, char **argv)
  {
    for (int i=1; i < argc; ++i) this->tokens.push_back(std::string(argv[i]));
  }

  const std::string& getCmdOption(const std::string &option) const
  {
    std::vector<std::string>::const_iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
      return *itr;
    }
    static const std::string empty_string("");
    return empty_string;
  }

  bool cmdOptionExists(const std::string &option) const{
    return std::find(this->tokens.begin(), this->tokens.end(), option)
      != this->tokens.end();
  }
private:
  std::vector <std::string> tokens;
};

// -------------------------------------------------------------
//  class App_Environment
// -------------------------------------------------------------
  class App_Environment
{
public:

  /// Default constructor.
  /** 
   * 
   * 
   * @param argc        number of program command line arguments
   * @param argv        command line arguments
   * @param help        help string to be displayed for the application
   * @param config_file name of the config file
   * 
   * @return 
   */
  App_Environment(int& argc, char **argv,
              const char* help,
	      const char* config_file);

  /// Destructor
  ~App_Environment(void);
private:
  // Command line parser
  CommandLineParser clparser;

  // environment for initializing the parallel libraries
  gridpack::parallel::Environment parenv;

  // Configuration file
  char config_file[100];
};



} // namespace app_environment
} // namespace gridpack

#endif
