/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   serial_io.hpp
 * @author Bruce Palmer
 * @date   2013-07-18 d3g293
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _serial_io_h_
#define _serial_io_h_

#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include "gridpack/parallel/distributed.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"

namespace gridpack {
namespace serial_io {

// -------------------------------------------------------------
// A set of classes to support output of information from buses
// and branches to standard output. Each bus or branch is
// responsible for creating a string that can be written to
// standard out. These modules then organize these sequentially
// and write them from process 0
// -------------------------------------------------------------

template <class _network>
class SerialBusIO {
  public:

  /**
   * Simple constructor
   * @param max_str_len the maximum string length written out by any bus
   * @param network the network for which output is desired
   */
  SerialBusIO(int max_str_len,
              boost::shared_ptr<_network> network)
  {
    p_GAgrp = network->communicator().getGroup();
    p_GA_type = NGA_Register_type(max_str_len);
    p_network = network;
    p_size = max_str_len;
  
    // Find total number of buses in network and create GAs for moving strings
    // around
    int nbus = p_network->totalBuses();
    int one = 1;
    p_stringGA = GA_Create_handle();
    GA_Set_data(p_stringGA,one,&nbus,p_GA_type);
    GA_Set_pgroup(p_stringGA, p_GAgrp);
    GA_Allocate(p_stringGA);
    p_maskGA = GA_Create_handle();
    GA_Set_data(p_maskGA,one,&nbus,C_INT);
    GA_Set_pgroup(p_maskGA, p_GAgrp);
    GA_Allocate(p_maskGA);
    p_useFile = false;
  }

  /**
   * Simple Destructor
   */
  ~SerialBusIO(void)
  {
    NGA_Deregister_type(p_GA_type);
    GA_Destroy(p_stringGA);
    GA_Destroy(p_maskGA);
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_useFile && p_fout->is_open()) p_fout->close();
    }
  }

  /**
   * Redirect output to a file instead of standard out
   * @param filename name of file that output goes to
   */
  void open(const char *filename)
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_useFile && p_fout->is_open()) p_fout->close();
      if (!p_useFile) p_fout.reset(new std::ofstream);
      p_fout->open(filename);
      p_useFile = true;
    }
  }

  /**
   * return IO stream
   * @return IO stream to file
   */
  boost::shared_ptr<std::ofstream> getStream()
  {
    return p_fout;
  }

  /**
   * Set IO stream to point to existing file
   * @param IO stream to file
   */
  void setStream(boost::shared_ptr<std::ofstream> stream)
  {
    p_fout = stream;
    p_useFile = true;
  }

  /**
   * Close file and redirect output to standard out 
   */
  void close()
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_fout->is_open()) p_fout->close();
    }
  }

  /**
   * Write output from buses to standard out
   * @param signal an optional character string used to control contents of
   *                output
   */
  void write(const char *signal = NULL)
  {
    if (p_useFile) {
      write(*p_fout, signal);
    } else {
      write(std::cout, signal);
    }
  }

  /**
   * Write single string to standard output. This is used to write headers for a
   * data listing. It is mostly a convenience function so that users do not have
   * to identify the head node
   * @param str character string containing the header
   */
  void header(const char *str)
  {
    if (p_useFile) {
      header(*p_fout, str);
    } else {
      header(std::cout, str);
    }
  }

  protected:

  /**
   * Write output from buses to standard out
   * @param out stream object for output
   * @param signal an optional character string used to control contents of
   *                output
   */
  void write(std::ostream & out, const char *signal = NULL)
  {
    int nBus = p_network->numBuses();
    char string[p_size];
    int nwrites = 0;
    int i;
    int one = 1;

    // Count up total strings being written from this processor
    for (i=0; i<nBus; i++) {
      if (p_network->getActiveBus(i) &&
          p_network->getBus(i)->serialWrite(string,p_size,signal)) {
        nwrites++;
      }
    }


    // Set up buffers to scatter strings to global buffer
    int **index;
    index = new int*[nwrites];
    int ones[nwrites];
    char strbuf[nwrites*p_size];
    char *ptr = strbuf;
    nwrites = 0;
    for (i=0; i<nBus; i++) {
      if (p_network->getActiveBus(i) &&
          p_network->getBus(i)->serialWrite(ptr,p_size,signal)) {
        index[nwrites] = new int;
        *(index[nwrites]) = p_network->getGlobalBusIndex(i);
        ones[nwrites] = 1;
        nwrites++;
        ptr += p_size;
      }
    }

    // Scatter data to global buffer and set mask array
    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      NGA_Scatter(p_stringGA,strbuf,index,nwrites);
      NGA_Scatter(p_maskGA,ones,index,nwrites);
    }
    GA_Pgroup_sync(p_GAgrp);
    for (i=0; i<nwrites; i++) {
      delete index[i];
    }
    delete [] index;

    // String data is now stored on global array. Process 0 now retrieves data
    // from each successive processor and writes it to standard out  
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      int nprocs = GA_Pgroup_nnodes(p_GAgrp);
      int lo, hi;
      for (i=0; i<nprocs; i++) {
        NGA_Distribution(p_maskGA, i, &lo, &hi);
        int ld = hi - lo + 1;
        // Figure out how many strings are coming from process i
        int imask[ld];
        NGA_Get(p_maskGA,&lo,&hi,imask,&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char iobuf[p_size*nwrites];
          index = new int*[nwrites];
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = new int;
              *(index[nwrites]) = j + lo;
              nwrites++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,index,nwrites);
          ptr = iobuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              out << ptr;
              ptr += p_size;
              delete index[nwrites];
              nwrites++;
            }
          }
          delete [] index;
        }
      }
    }
    GA_Pgroup_sync(p_GAgrp);
  }

  /**
   * Write single string to standard output. This is used to write headers for a
   * data listing. It is mostly a convenience function so that users do not have
   * to identify the head node
   * @param out stream object for output
   * @param str character string containing the header
   */
  void header(std::ostream & out, const char *str) const
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      out << str;
    }
  }

  private:
    int p_GA_type;
    boost::shared_ptr<_network> p_network;
    int p_stringGA;
    int p_maskGA;
    int p_size;
    bool p_useFile;
    boost::shared_ptr<std::ofstream> p_fout;
    int p_GAgrp;
};

template <class _network>
class SerialBranchIO {
  public:
  /**
   * Simple constructor
   * @param max_str_len the maximum string length written out by any branch
   * @param network the network for which output is desired
   */
  SerialBranchIO(int max_str_len,
                 boost::shared_ptr<_network> network)
  {
    p_GAgrp = network->communicator().getGroup();
    p_GA_type = NGA_Register_type(max_str_len);
    p_network = network;
    p_size = max_str_len;

    // Find total number of branches in network and create GAs for moving strings
    // around
    int nbranch = p_network->totalBranches();
    int one = 1;
    p_stringGA = GA_Create_handle();
    GA_Set_data(p_stringGA,one,&nbranch,p_GA_type);
    GA_Set_pgroup(p_stringGA, p_GAgrp);
    GA_Allocate(p_stringGA);
    p_maskGA = GA_Create_handle();
    GA_Set_data(p_maskGA,one,&nbranch,C_INT);
    GA_Set_pgroup(p_maskGA, p_GAgrp);
    GA_Allocate(p_maskGA);
    p_useFile = false;
  }

  /**
   * Simple Destructor
   */
  ~SerialBranchIO(void)
  {
    NGA_Deregister_type(p_GA_type);
    GA_Destroy(p_stringGA);
    GA_Destroy(p_maskGA);
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_useFile && p_fout->is_open()) p_fout->close();
    }
  }

  /**
   * Redirect output to a file instead of standard out
   * @param filename name of file that output goes to
   */
  void open(const char *filename)
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_useFile && p_fout->is_open()) p_fout->close();
      if (!p_useFile) p_fout.reset(new std::ofstream);
      p_fout->open(filename);
      p_useFile = true;
    }
  }

  /**
   * return IO stream
   * @return IO stream to file
   */
  boost::shared_ptr<std::ofstream> getStream()
  {
    return p_fout;
  }

  /**
   * Set IO stream to point to existing file
   * @param IO stream to file
   */
  void setStream(boost::shared_ptr<std::ofstream> stream)
  {
    p_fout = stream;
    p_useFile = true;
  }

  /**
   * Close file and redirect output to standard out 
   */
  void close()
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_useFile && p_fout->is_open()) p_fout->close();
    }
  }

  /**
   * Write output from branches to standard out
   * @param signal an optional character string used to control contents of
   *                output
   */
  void write(const char *signal = NULL)
  {
    if (p_useFile) {
      write(*p_fout, signal);
    } else {
      write(std::cout, signal);
    }
  }

  /**
   * Write single string to standard output. This is used to write headers for a
   * data listing. It is mostly a convenience function so that users do not have
   * to identify the head node
   * @param str character string containing the header
   */
  void header(const char *str)
  {
    if (p_useFile) {
      header(*p_fout, str);
    } else {
      header(std::cout, str);
    }
  }

  protected:

  /**
   * Write output from branches to standard out
   * @param out stream object for output
   * @param signal an optional character string used to control contents of
   *                output
   */
  void write(std::ostream & out, const char *signal = NULL)
  {
    int nBranch = p_network->numBranches();
    char string[p_size];
    int nwrites = 0;
    int i;
    int one = 1;

    // Count up total strings being written from this processor
    for (i=0; i<nBranch; i++) {
      if (p_network->getActiveBranch(i) &&
          p_network->getBranch(i)->serialWrite(string,p_size,signal)) nwrites++;
    }

    // Set up buffers to scatter strings to global buffer
    int **index;
    index = new int*[nwrites];
    int ones[nwrites];
    char strbuf[nwrites*p_size];
    char *ptr = strbuf;
    nwrites = 0;
    for (i=0; i<nBranch; i++) {
      if (p_network->getActiveBranch(i) &&
          p_network->getBranch(i)->serialWrite(ptr,p_size,signal)) {
        index[nwrites] = new int;
        *(index[nwrites]) = p_network->getGlobalBranchIndex(i);
        ones[nwrites] = 1;
        nwrites++;
        ptr += p_size;
      }
    }

    // Scatter data to global buffer and set mask array
    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      NGA_Scatter(p_stringGA,strbuf,index,nwrites);
      NGA_Scatter(p_maskGA,ones,index,nwrites);
    }
    GA_Pgroup_sync(p_GAgrp);
    for (i=0; i<nwrites; i++) {
      delete index[i];
    }
    delete [] index;

    // String data is now stored on global array. Process 0 now retrieves data
    // from each successive processor and writes it to standard out  
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      int nprocs = GA_Pgroup_nnodes(p_GAgrp);
      int lo, hi;
      for (i=0; i<nprocs; i++) {
        NGA_Distribution(p_maskGA, i, &lo, &hi);
        int ld = hi - lo + 1;
        // Figure out how many strings are coming from process i
        int imask[ld];
        NGA_Get(p_maskGA,&lo,&hi,imask,&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char iobuf[p_size*nwrites];
          index = new int*[nwrites];
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = new int;
              *(index[nwrites]) = j + lo;
              nwrites++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,index,nwrites);
          ptr = iobuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              out << ptr;
              ptr += p_size;
              delete index[nwrites];
              nwrites++;
            }
          }
          delete [] index;
        }
      }
    }
    GA_Pgroup_sync(p_GAgrp);
  }

  /**
   * Write single string to standard output. This is used to write headers for a
   * data listing. It is mostly a convenience function so that users do not have
   * to identify the head node
   * @param out stream object for output
   * @param str character string containing the header
   */
  void header(std::ostream & out, const char *str) const
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      out << str;
    }
  }

  private:
    int p_GA_type;
    boost::shared_ptr<_network> p_network;
    int p_stringGA;
    int p_maskGA;
    int p_size;
    bool p_useFile;
    boost::shared_ptr<std::ofstream> p_fout;
    int p_GAgrp;
};

}   // serial_io
}   // gridpack
#endif  // _serial_io_h_
