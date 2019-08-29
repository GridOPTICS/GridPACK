// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   bus_table.hpp
 * @author Bruce Palmer
 * @date   2016-12-20 09:25:03 d3g293
 * 
 * @brief  
 * This is a utility that is designed to parse an external file representing a
 * collection of data that maps to individual buses and make that data available
 * on the processors containing the buses. The file has the logical structure
 * busID tag value0 value1 .... valueN-1
 * for multiple buses. BusID is an integer corresponding to the original bus
 * index, tag is a 1 or 2 character string that maps to some object on the bus
 * and value0, value1, etc. are different values that map to each bus object.
 * Each bus object has the same number of values mapping to it.
 * 
 */

// -------------------------------------------------------------

#ifndef _bus_table_hpp_
#define _bus_table_hpp_

#define NUM_SLICES 10

#include <ga.h>
#include <boost/unordered_map.hpp>
#include "gridpack/parser/hash_distr.hpp"
#include "gridpack/utilities/exception.hpp"
#include "gridpack/utilities/string_utils.hpp"

namespace gridpack {
namespace bus_table {

typedef struct {
  int order;
  char tag[3];
} table_t;

// -------------------------------------------------------------
//  class HashDistribution
// -------------------------------------------------------------
template <typename _network>
 class BusTable {

public:
  typedef _network NetworkType;
  typedef boost::shared_ptr<NetworkType> NetworkPtr;

  // Default constructor
  // @param network network over which hash distribution function extends
  BusTable(const boost::shared_ptr<_network> network)
    : p_network(network)
  {
    p_comm = p_network->communicator();
    p_has_data = false;
  }

  // Default destructor
  ~BusTable(void)
  {
    if (p_has_data) GA_Destroy(p_GA_data);
  }

  // Store contents of file in a distributed bus table
  bool readTable(std::string &filename)
  {
    if (p_has_data) GA_Destroy(p_GA_data);
    int me = p_comm.rank();
    MPI_Comm mpi_comm = static_cast<MPI_Comm>(p_comm);
    p_GAgrp = p_comm.getGroup();
    int ifound = 1;
    int nval = 0;
    int nline = 0;
    gridpack::utility::StringUtils util;
    p_headers.clear();
    if (me == 0) {
      // Cheesy hack to find out how many lines are in the file. Just open
      // the file and read all lines, then close it and open it again.
      std::ifstream input;
      input.open(filename.c_str());
      if (!input.is_open()) {
        std::cout<<"File "<<filename<<" not found by readTable"<<std::endl;
        ifound = 0;
      } else {
        std::string line;
        // Check to see if any lines are comments
        if (!std::getline(input,line).eof()) {
          util.trim(line);
          bool found = true;
          while (line[0] == '#') {
            // Line is a comment. Store as a string for use elsewhere after
            // removing the comment character (#)
            line[0] = ' ';
            util.trim(line);
            p_headers.push_back(line);
            found = !std::getline(input,line).eof();
          }
          if (found) {
            util.trim(line);
            std::vector<std::string> split_line;
            boost::split(split_line,line,boost::algorithm::is_any_of(" "),
                boost::token_compress_on);
            nval = split_line.size() - 2;
            if (nval > 0) {
              nline = 1;
              while(!std::getline(input,line).eof()) {
                nline++;
              }
            } else {
              ifound = 0;
            }
          }
        } else {
          ifound = 0;
        }
      }
      input.close();
      // At this point, process 0 knows whether any data has been found,
      // how many values there are per line, and how many lines there are.
      // Based on this information, it is possible to set up some distributed
      // data structures to hold all this information. First, broadcast these
      // values to all processors.
      int irecv;
      MPI_Bcast(&ifound,1,MPI_INT,0,mpi_comm);
      if (ifound == 1) {
        MPI_Bcast(&nval,1,MPI_INT,0,mpi_comm);
        MPI_Bcast(&nline,1,MPI_INT,0,mpi_comm);
      } else {
        return false;
      }
    } else {
      MPI_Bcast(&ifound,1,MPI_INT,0,mpi_comm);
      if (ifound == 1) {
        MPI_Bcast(&nval,1,MPI_INT,0,mpi_comm);
        MPI_Bcast(&nline,1,MPI_INT,0,mpi_comm);
      } else {
        return false;
      }
    }
    p_nobjs = nline;
    p_nvals = nval;
    
    // File has been found and data is present. Create a global array to hold
    // the data and reread the file, this time parsing all data.
    // Also broadcast and header data that might be available.
    if (nval > 0 && nline > 0) {
      // Distribute headers
      int nheaders = 0;
      if (me == 0) {
        nheaders = p_headers.size();
      }
      MPI_Bcast(&nheaders,1,MPI_INT,0,mpi_comm);
      int *hsize;
      hsize = new int[nheaders];
      int i;
      if (me == 0) {
        for (i=0; i<nheaders; i++) {
          hsize[i] = p_headers[i].length();
        }
      }
      MPI_Bcast(hsize,nheaders,MPI_INT,0,mpi_comm);
      for (i=0; i<nheaders; i++) {
        char *str;
        str = new char[hsize[i]];
        int j;
        if (me == 0) {
          for (j=0; j<hsize[i]; j++) {
            str[j] = (p_headers[i])[j];
          }
        }
        MPI_Bcast(str,hsize[i],MPI_CHAR,0,mpi_comm);
        if (me != 0) {
          p_headers.push_back(str);
        }
        delete [] str;
      }
      delete [] hsize;
      gridpack::utility::StringUtils util;
      // Create GA to hold all data
      p_GA_data = GA_Create_handle();
      int ndim = 1;
      int dims = nline*nval;
      GA_Set_data(p_GA_data,ndim,&dims,C_DBL);
      GA_Set_pgroup(p_GA_data,p_GAgrp);
      GA_Allocate(p_GA_data);
      std::vector<table_t> order;
      std::vector<int> bus_id;
      if (me == 0) {
        int i;
        int ncnt = 0;
        int lo, hi;
        lo = 0;
        // allocate local buffer for values
        double *buf = new double[NUM_SLICES*nval];
        double *ptr = buf;
        std::ifstream input;
        input.open(filename.c_str());
        std::vector<std::string> split_line;
        std::string line;
        while(!std::getline(input,line).eof()) {
          while (line[0] == '#') {
            // Line is a comment. Store as a string for use elsewhere after
            // removing the comment character (#)
            std::getline(input,line);
          }
          util.trim(line);
          boost::split(split_line,line,boost::algorithm::is_any_of(" "),
              boost::token_compress_on);
          table_t data;
          data.order = ncnt;
          std::string tag = split_line[1];
          std::string new_tag = util.clean2Char(tag);
          strncpy(data.tag,new_tag.c_str(),2);
          data.tag[2] = '\0';
          bus_id.push_back(atoi(split_line[0].c_str()));
          order.push_back(data);
          for (i=0; i<nval; i++) {
            ptr[i] = atof(split_line[i+2].c_str());
          }
          ptr += nval;
          ncnt++;
          if (ncnt%NUM_SLICES == 0) {
            hi = ncnt*nval-1;
            NGA_Put(p_GA_data,&lo,&hi,buf,&nval);
            ptr = buf;
            lo = ncnt*nval;
          }
        }
        if (ncnt != nline) {
          printf("Mismatch parsing data in table %s\n",filename.c_str());
        }
        if (ncnt%NUM_SLICES != 0) {
          hi = nline*nval-1;
          NGA_Put(p_GA_data,&lo,&hi,buf,&nval);
        }
        input.close();
        delete [] buf;
      }
      // Synchronize across all processors so data is in a consistent state
      p_comm.sync();
      // Create hash table to match data with processors holding corresponding
      // buses
      gridpack::hash_distr::HashDistribution<NetworkType,table_t,table_t>
        *hash;
      hash = new gridpack::hash_distr::HashDistribution
        <NetworkType,table_t,table_t>(p_network);
      hash->distributeBusValues(bus_id,order);
      delete hash;
      p_local_idx.clear();
      p_tags.clear();
      p_order.clear();
      for (i=0; i<bus_id.size(); i++) {
        p_local_idx.push_back(bus_id[i]);
        std::string tmp(order[i].tag);
        std::string new_tag = util.clean2Char(tmp);
        p_tags.push_back(new_tag);
        p_order.push_back(order[i].order);
      }
    }
    p_has_data = true;
    return true;
  }

  /**
   * Return a list of local indices of all buses in the table
   * @param indices array holding local indices
   */
  void getLocalIndices(std::vector<int> &indices) {
    indices.clear();
    int nval = p_local_idx.size();
    int i;
    for (i=0; i<nval; i++) {
      indices.push_back(p_local_idx[i]);
    }
  }

  /**
   * Return a list of tags for all objects on buses in the table
   * @param tags array holding tags
   */
  void getTags(std::vector<std::string> &tags) {
    tags.clear();
    int nval = p_tags.size();
    int i;
    for (i=0; i<nval; i++) {
      tags.push_back(p_tags[i]);
    }
  }

  /**
   * Return an array of values for the ith column of data in the table
   * The vector of values can be matched to the corresponding bus object
   * using the array returned by getLocalIndices and getTags
   * @param idx index of column of data
   * @param values vector containing return values.
   */
  void getValues(int idx, std::vector<double> &values) {
    values.clear();
    if (idx < 0 || idx >= p_nvals) {
      printf("Requested data is out of range of bus table (%d not in [0,%d])\n",
          idx,p_nvals-1);
    }
    //construct array of indices to retrieve data
    int nsize = p_local_idx.size();
    std::vector<double> v(nsize);
    std::vector<int> array(nsize);
    std::vector<int*> subscript(nsize);
    int *ptr = &array[0];
    int i;
    for (i=0; i<nsize; i++) {
      array[i] = p_nvals*p_order[i] + idx;
      subscript[i] = ptr;
      ptr++;
    }
    NGA_Gather(p_GA_data, &v[0], &subscript[0], nsize);
    p_comm.sync();
    for (i=0; i<nsize; i++) {
      values.push_back(v[i]);
    }
  }

  void getHeaders(std::vector<std::string> &headers) {
    headers.clear();
    int i;
    for (i=0; i<p_headers.size(); i++) {
      headers.push_back(p_headers[i]);
    }
  }

  /**
   * Return the number of columns in the table
   */
  int getNumColumns()
  {
    return p_nvals;
  }
private:

  NetworkPtr p_network;
  
  gridpack::parallel::Communicator p_comm;
  int p_GA_data;
  bool p_has_data;
  int p_nvals;
  int p_nobjs;

  std::vector<int> p_local_idx;
  std::vector<std::string> p_tags;
  std::vector<int> p_order;
  std::vector<std::string> p_headers;

  int p_GAgrp;
};


} // namespace bus_table
} // namespace gridpack

#endif

