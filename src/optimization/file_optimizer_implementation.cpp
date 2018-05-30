// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   file_optimizer_implementation.cpp
 * @author William A. Perkins
 * @date   2016-12-08 14:59:36 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December  7, 2016 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include "file_optimizer_implementation.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class FileOptimizerImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// FileOptimizerImplementation:: constructors / destructor
// -------------------------------------------------------------
FileOptimizerImplementation::FileOptimizerImplementation(const parallel::Communicator& comm)
  : OptimizerImplementation(comm)
{
  
}

FileOptimizerImplementation::~FileOptimizerImplementation(void)
{
}

// -------------------------------------------------------------
// FileOptimizerImplementation::p_configure
// -------------------------------------------------------------
void 
FileOptimizerImplementation::p_configure(utility::Configuration::CursorPtr props)
{
  p_outputName = props->get("File", p_outputName);
  if (!p_outputName.empty()) {
    int me(this->processor_rank());
    p_outputName += boost::str(boost::format("%04d") % me);
  } else {
    p_outputName = p_temporaryFileName();
  }
  p_runMaybe = props->get("Run", true);
}

// -------------------------------------------------------------
// FileOptimizerImplementation::p_setFilename
// -------------------------------------------------------------
void
FileOptimizerImplementation::p_setFilename(std::string file)
{
  p_outputName = file;
}

// -------------------------------------------------------------
// FileOptimizerImplementation::p_temporaryFile
// -------------------------------------------------------------
std::string
FileOptimizerImplementation::p_temporaryFileName(void)
{
  using namespace boost::filesystem;
  std::string pattern = 
    boost::str(boost::format("%s_%d") % "gridpack%%%%" % this->processor_rank());
  path model(pattern);
  path tmp(temp_directory_path());
  tmp /= unique_path(model);

  boost::system::error_code ec;
  file_status istat = status(tmp);
  std::string result(tmp.c_str());
  return result;
  
}

// -------------------------------------------------------------
// FileOptimizerImplementation::p_solve
// -------------------------------------------------------------
void
FileOptimizerImplementation::p_solve(const p_optimizeMethod& method)
{
  std::string tmpname(p_outputName);
  std::ofstream tmp;
  tmp.open(tmpname.c_str());
  if (!tmp) {
    std::string msg("Cannot open temporary file: ");
    msg += tmpname.c_str();
    throw gridpack::Exception(msg);
  }
  p_write(method, tmp);
  tmp.close();
}

// -------------------------------------------------------------
// FileOptimizerImplementation
// -------------------------------------------------------------



} // namespace optimization
} // namespace gridpack
