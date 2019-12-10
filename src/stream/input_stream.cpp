/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

#include <iostream>
#include "gridpack/stream/input_stream.hpp"

/**
 * Constructor
 */
gridpack::stream::InputStream::InputStream()
{
  p_srcFile = false;
  p_srcChannel = false;
  p_isOpen = false;
}

/**
 * Destructor
 */
gridpack::stream::InputStream::~InputStream()
{
  if (p_isOpen) {
    if (p_srcFile) {
      p_fout.close();
    } else if (p_srcChannel) {
#ifdef USE_GOSS
#endif
    }
  }
}

/**
 * Open a file
 * @param filename name of file to open
 * @return true of file is successfully opened
 */
bool gridpack::stream::InputStream::openFile(std::string file)
{
  bool ret = false;
  p_fout.open(file.c_str());
  if (p_fout.is_open()) {
    p_srcFile = true;
    p_isOpen = true;
    ret = true;
  }
  return ret;
}

#ifdef USE_GOSS
/**
 * Open GOSS channel
 * @param topic topic that that is used by server
 * @param URI channel URI
 * @param username server username
 * @param passwd server password
 */
bool gridpack::stream::InputStream::openChannel(const char *topic,
    const char *URI, const char *username, const char *passwd)
{
  bool ret = false;
  return ret;
}
#endif

void gridpack::stream::InputStream::close()
{
  if (p_isOpen) {
    if (p_srcFile) {
      p_fout.close();
    } else if (p_srcChannel) {
#ifdef USE_GOSS
#endif
    }
  } else {
    std::cout<<"No input stream open when calling close"<<std::endl;
  }
}

/**
 * Get next line from streaming source
 * @param line next line from streaming source
 * @return true if another line is found
 */
bool gridpack::stream::InputStream::nextLine(std::string &line)
{
  bool ret = false;
  if (p_isOpen) {
    if (p_srcFile) {
      ret = std::getline(p_fout, line).good();
    } else {
#ifdef USE_GOSS
#endif
    }
  } else {
  }
  return ret;
}

/**
 * Report if a stream is currently open
 * @return true if stream is open
 */
bool gridpack::stream::InputStream::isOpen()
{
  return p_isOpen;
}
