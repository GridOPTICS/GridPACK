/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#ifndef _input_stream_h
#define _input_stream_h

#include <string>
#include <vector>
#include <fstream>

// Simple class to open a stream of data and read lines

namespace gridpack{
namespace stream{

class InputStream {
public:

  /**
   * Constructor
   */
  InputStream();

  /**
   * Destructor
   */
  ~InputStream();

  /**
   * Open a file
   * @param filename name of file to open
   * @return true of file is successfully opened
   */
  bool openFile(std::string file);

  /**
   * Parse a vector of strings representing a file
   * @param fileVec a vector of strings containing lines in a file
   */
  bool openStringVector(const std::vector<std::string> &fileVec);

  /**
   * Close a file or other input stream
   */
  void close();

  /**
   * Get next line from streaming source
   * @param line next line from streaming source
   * @return true if another line is found
   */
  bool nextLine(std::string &line);
  
  /**
   * Report if a stream is currently open
   * @return true if stream is open
   */
  bool isOpen();

private:

  std::ifstream p_fout;

  bool p_srcFile;

  bool p_srcVector;

  bool p_isOpen;

  std::vector<std::string> p_fileVector;

  std::vector<std::string>::iterator p_fileIterator;
};


}    // stream
}    // gridpack


#endif // _input_stream_h
