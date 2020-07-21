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

#include <boost/mpi.hpp>
#include "gridpack/utilities/uncopyable.hpp"


namespace gridpack {

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
//  class Environment
// -------------------------------------------------------------
  class Environment : private utility::Uncopyable
{
public:

  /// Default constructor.
  /** 
   * 
   * 
   * @param argc        number of program command line arguments
   * @param argv        command line arguments
   * @param help        help string to be displayed for the application
   * 
   * @return 
   */
  Environment(int argc, char **argv,
              const char* help);

  Environment(int argc, char **argv,
              const char* help,
	      const long int& ma_stack,
              const long int& ma_heap);

  Environment(int argc, char **argv);


  /// Destructor
  ~Environment(void);

  /**
   * return next token after option
   */
  const std::string getCmdOption(std::string &option);

protected:
  boost::mpi::environment p_boostEnv;

private:
  // Command line parser
  CommandLineParser clparser;

  // Stack and heap for GA
  long int pma_stack;
  long int pma_heap;

  void PrintHelp(char **argv,const char* help);

  void PrintStatus();
};

} // namespace gridpack

#endif
