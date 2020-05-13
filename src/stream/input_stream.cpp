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
  p_srcVector = false;
  p_isOpen = false;
}

/**
 * Destructor
 */
gridpack::stream::InputStream::~InputStream()
{
  if (p_isOpen) {
    p_fout.close();
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

/**
 * Parse a vector of strings representing a file
 * @param fileVec a vector of strings containing lines in a file
 */
bool gridpack::stream::InputStream::openStringVector(const std::vector<std::string> &fileVec)
{
  if (fileVec.size() == 0) return false;
  p_fileVector = fileVec;
  p_srcVector = true;
  p_isOpen = true;
  return true;
}

void gridpack::stream::InputStream::close()
{
  if (p_isOpen) {
    if (p_srcFile) {
      p_fout.close();
      p_srcFile = false;
    } else if (p_srcVector) {
      p_fileVector.clear();
      p_srcVector = false;
    }
    p_isOpen = false;
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
      if (p_fileIterator != p_fileVector.end()) {
        line = *p_fileIterator;
        ret = true;
        p_fileIterator++;
      } else {
        line.clear();
      }
    }
  } else {
    std::cout<<"No input stream is open when calling nextLine"<<std::endl;
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
