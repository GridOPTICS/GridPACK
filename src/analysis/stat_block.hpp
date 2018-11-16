// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   stat_block.hpp
 * @author Bruce Palmer
 * @date   2018-09-11
 * 
 * @brief  
 * This is a utility that is designed to allow users to create a large
 * distributed table of data that can subsequently be use for statistical
 * analysis. Values in the table are masked so that only values that have been
 * deemed relevant according to some criteria are included in the analysis.
 * 
 */

// -------------------------------------------------------------

#ifndef _stat_block_hpp_
#define _stat_block_hpp_

#include <ga.h>
#include <map>
#include <vector>
#include "gridpack/parallel/communicator.hpp"

namespace gridpack {
namespace analysis {

class StatBlock {
private:
  // Simple data struct to keep track of row labels
  typedef struct {
                 int gidx;
                 int idx1;
                 int idx2;
                 char tag[3];
  } index_set;

public:
  /**
   * Constructor
   * @param comm communicator on which StatBlock is defined
   * @param nrows number of rows in data array
   * @param ncols number of columns in data array
   */
  StatBlock(const parallel::Communicator &comm, int nrows, int ncols);

  /**
   * Default destructor
   */
  ~StatBlock(void);

  /**
   * Add a column of data to the stat block
   * @param idx index of column
   * @param vals vector of column values
   * @param mask vector of mask values
   */
  void addColumnValues(int idx, std::vector<double> vals, std::vector<int> mask);

  /**
   * Add index and device tag that can be used to label rows
   * @param indices vector of indices
   * @param tags  vector of character tags
   */
  void addRowLabels(std::vector<int> indices, std::vector<std::string> tags);

  /**
   * Add two branch indices and device tag that can be used to label rows
   * @param idx1 vector of index 1
   * @param idx2 vector of index 2
   * @param tags  vector of character tags
   */
  void addRowLabels(std::vector<int> idx1, std::vector<int> idx2,
      std::vector<std::string> tags);

  /**
   * Add the minimum allowed value per row
   * @param max vector containing minimum value for each row
   */
  void addRowMinValue(std::vector<double> min);

  /**
   * Add the maximum allowed value per row
   * @param max vector containing maximum value for each row
   */
  void addRowMaxValue(std::vector<double> max);

  /**
   * Write out file containing mean value and RMS deviation for values in table
   * @param filename name of file containing results
   * @param mval only include values with this mask value
   * @param flag if false, do not include tag ids in output
   */
  void writeMeanAndRMS(std::string filename, int mval=1, bool flag = true);

  /**
   * Write out file containing Min an Max values in table for each row
   * @param filename name of file containing results
   * @param mval only include values with this mask value or greater
   * @param flag if false, do not include tag ids in output
   */
  void writeMinAndMax(std::string filename, int mval=1, bool flag = true);

  /**
   * Write out file containing number of mask entries at each row that
   * correspond to a given value
   * @param filename name of file containing results
   * @param mval count number of times this mask value occurs
   * @param flag if false, do not include tag ids in output
   */
  void writeMaskValueCount(std::string filename, int mval, bool flag = true);

  /**
   * Sum up the values in the columns and print the result as a function
   * of column index
   * @param filename name of file containing results
   * @param mval only include values with this mask value or greater
   */
  void sumColumnValues(std::string filename, int mval=1);
private:

  int p_data;
  int p_mask;
  int p_type;
  int p_tags;
  int p_bounds;

  bool p_max_bound;
  bool p_min_bound;
  bool p_branch_flag;
   
  int p_nrows;
  int p_ncols;

  int p_nprocs;
  int p_me;
  int p_GAgrp;

  MPI_Comm p_comm;

};


} // namespace gridpack
} // namespace analysis

#endif
