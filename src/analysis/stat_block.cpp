// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   stat_block.cpp
 * @author Bruce Palmer
 * @date   2019-11-14 07:11:55 d3g096
 * 
 * @brief  
 * This is a utility that is designed to allow users to create a large
 * distributed table of data that can subsequently be use for statistical
 * analysis. Values in the table are masked so that only values that have been
 * deemed relevant according to some criteria are included in the analysis.
 * 
 */

// -------------------------------------------------------------

#include "gridpack/analysis/stat_block.hpp"
#include "gridpack/utilities/string_utils.hpp"

#define stb gridpack::analysis::StatBlock

#define BLOCKSIZE 100

#include <fstream>

/**
 * Constructor
 * @param comm communicator on which StatBlock is defined
 * @param nrows number of rows in data array
 * @param ncols number of columns in data array
 */
stb::StatBlock(const parallel::Communicator &comm, int nrows, int ncols)
{
  int one = 1;
  int two = 2;
  int dims[2];
  int chunk[2];
  p_min_bound = false;
  p_max_bound = false;
  p_nrows = nrows; 
  p_ncols = ncols; 
  p_nprocs = comm.size();
  p_me = comm.rank();
  p_comm = static_cast<MPI_Comm>(comm);
  p_GAgrp = comm.getGroup();
  p_branch_flag = false;

  // Create data and mask arrays
  dims[0] = nrows;
  dims[1] = ncols;
  chunk[0] = -1;
  chunk[1] = -1;

  p_data = GA_Create_handle();
  GA_Set_data(p_data,two,dims,C_DBL);
  GA_Set_chunk(p_data,chunk);
  GA_Set_pgroup(p_data,p_GAgrp);
  GA_Allocate(p_data);

  p_mask = GA_Create_handle();
  GA_Set_data(p_mask,two,dims,C_INT);
  GA_Set_chunk(p_mask,chunk);
  GA_Set_pgroup(p_mask,p_GAgrp);
  GA_Allocate(p_mask);


  p_type = NGA_Register_type(sizeof(index_set));
  p_tags = GA_Create_handle();
  GA_Set_data(p_tags,one,&p_nrows,p_type);
  GA_Set_pgroup(p_tags,p_GAgrp);
  GA_Allocate(p_tags);

  dims[0] = p_nrows;
  dims[1] = 2;

  p_bounds = GA_Create_handle();
  GA_Set_data(p_bounds,two,dims,C_DBL);
  GA_Set_pgroup(p_bounds,p_GAgrp);
  GA_Allocate(p_bounds);
}

/**
 * Default destructor
 */
stb::~StatBlock(void)
{
  NGA_Deregister_type(p_type);
  GA_Destroy(p_data);
  GA_Destroy(p_mask);
  GA_Destroy(p_tags);
  GA_Destroy(p_bounds);
}

/**
 * Add a column of data to the stat block
 * @param idx index of column
 * @param vals vector of column values
 * @param mask vector of mask values
 */
void stb::addColumnValues(int idx, std::vector<double> vals, std::vector<int> mask)
{
  if (idx <p_ncols && idx >= 0) {
    int lo[2];
    int hi[2];
    int ld = 1;
    lo[0] = 0;
    hi[0] = p_nrows-1;
    lo[1] = idx;
    hi[1] = idx;
    NGA_Put(p_data,lo,hi,&vals[0],&ld);
    NGA_Put(p_mask,lo,hi,&mask[0],&ld);
  } else {
    printf("IDX: %d NCOLS: %d\n",idx,p_ncols);
    // TODO: Some kind of error
  }
}

/**
 * Add index and device tag that can be used to label rows
 * @param indices vector of indices
 * @param tags  vector of character tags
 */
void stb::addRowLabels(std::vector<int> indices, std::vector<std::string> tags)
{
  std::vector<stb::index_set> tagvec;
  gridpack::utility::StringUtils util;
  if (indices.size() != p_nrows) {
    printf("NROWS: %d indices.size: %d\n",p_nrows,(int)indices.size());
    // TODO: Some kind of error
  }
  if (p_nrows != tags.size()) {
    // TODO: Some kind of error
    printf("NROWS: %d tags.size: %d\n",p_nrows,(int)tags.size());
  }
  int i;
  for (i=0; i<indices.size(); i++) {
    stb::index_set tag;
    tag.gidx = i+1;
    tag.idx1 = indices[i];
    std::string ctk;
    if (tags[i].size() > 1) {
      ctk = tags[i].substr(0,2);
    } else {
      ctk = tags[i];
    }
    util.clean2Char(ctk);
    strncpy(tag.tag,ctk.c_str(),2);
    tag.tag[2] = '\0';
    tagvec.push_back(tag);
  }
  int lo = 0;
  int hi = p_nrows-1;
  int one = 1;
  NGA_Put(p_tags,&lo,&hi,&tagvec[0],&one);
}

/**
 * Add two branch indices and device tag that can be used to label rows
 * @param idx1 vector of index 1
 * @param idx2 vector of index 2
 * @param tags  vector of character tags
 */
void stb::addRowLabels(std::vector<int> idx1, std::vector<int> idx2,
    std::vector<std::string> tags)
{
  p_branch_flag = true;
  std::vector<stb::index_set> tagvec;
  gridpack::utility::StringUtils util;
  if (idx1.size() != p_nrows) {
    printf("NROWS: %d idx1.size: %d\n",p_nrows,(int)idx1.size());
    // TODO: Some kind of error
  }
  if (idx2.size() != p_nrows) {
    printf("NROWS: %d idx2.size: %d\n",p_nrows,(int)idx2.size());
    // TODO: Some kind of error
  }
  if (p_nrows != tags.size()) {
    printf("NROWS: %d tags.size: %d\n",p_nrows,(int)tags.size());
    // TODO: Some kind of error
  }
  int i;
  for (i=0; i<idx1.size(); i++) {
    stb::index_set tag;
    tag.gidx = i+1;
    tag.idx1 = idx1[i];
    tag.idx2 = idx2[i];
    std::string ctk;
    if (tags[i].size() > 1) {
      ctk = tags[i].substr(0,2);
    } else {
      ctk = tags[i];
    }
    util.clean2Char(ctk);
    strncpy(tag.tag,ctk.c_str(),2);
    tag.tag[2] = '\0';
    tagvec.push_back(tag);
  }
  if (p_nrows != tags.size()) {
    // TODO: Some kind of error
  }
  for (i=0; i<idx1.size(); i++) {
    stb::index_set tag;
    tag.gidx = i+1;
    tag.idx1 = idx1[i];
    tag.idx2 = idx2[i];
    std::string ctk;
    if (tags[i].size() > 1) {
      ctk = tags[i].substr(0,2);
    } else {
      ctk = tags[i];
    }
    util.clean2Char(ctk);
    strncpy(tag.tag,ctk.c_str(),2);
    tag.tag[2] = '\0';
    tagvec.push_back(tag);
  }
  int lo = 0;
  int hi = p_nrows-1;
  int one = 1;
  NGA_Put(p_tags,&lo,&hi,&tagvec[0],&one);
}

/**
 * Add the minimum allowed value per row
 * @param max vector containing minimum value for each row
 */
void stb::addRowMinValue(std::vector<double> min)
{
  if (min.size() != p_nrows) {
    printf("NROWS: %d min.size: %d\n",p_nrows,(int)min.size());
    // TODO: Some kind of error
  }
  int one = 1;
  int lo[2], hi[2];
  lo[0] = 0;
  hi[0] = p_nrows-1;
  lo[1] = 0;
  hi[1] = 0;
  NGA_Put(p_bounds,lo,hi,&min[0],&one);
  p_min_bound = true;
}

/**
 * Add the maximum allowed value per row
 * @param max vector containing maximum value for each row
 */
void stb::addRowMaxValue(std::vector<double> max)
{
  if (max.size() != p_nrows) {
    printf("NROWS: %d min.size: %d\n",p_nrows,(int)max.size());
    // TODO: Some kind of error
  }
  int one = 1;
  int lo[2], hi[2];
  lo[0] = 0;
  hi[0] = p_nrows-1;
  lo[1] = 1;
  hi[1] = 1;
  NGA_Put(p_bounds,lo,hi,&max[0],&one);
  p_max_bound = true;
}

/**
 * Write out file containing mean value and RMS deviation for values in table
 * @param filename name of file containing results
 * @param mval only include values with this mask value or greater
 * @param flag if false, do not include tag ids in output
 */
void stb::writeMeanAndRMS(std::string filename, int mval, bool flag)
{
  GA_Pgroup_sync(p_GAgrp);
  int zero = 0;
  int one = 1;
  int two = 2;
  int g_cnt = GA_Create_handle();
  NGA_Set_data(g_cnt,one,&one,C_INT);
  GA_Set_pgroup(g_cnt,p_GAgrp);
  NGA_Allocate(g_cnt);
  GA_Zero(g_cnt);
  int nblock = p_nrows/BLOCKSIZE;
  while (p_nrows > nblock*BLOCKSIZE) {
    nblock++;
  }
  int itask = NGA_Read_inc(g_cnt,&zero,(long)one);
  int dims[2];
  int g_buf = GA_Create_handle();
  dims[0] = p_nrows;
  dims[1] = 3;
  NGA_Set_data(g_buf,two,dims,C_DBL);
  GA_Set_pgroup(g_buf,p_GAgrp);
  NGA_Allocate(g_buf);
  std::vector<double> vavg, vavg2, vdiff2;
  int lo[2];
  int hi[2];
  int ld;
  int i, j; 
  while (itask < nblock) {
    // Buffers to hold blocks of data
    int blocksize;
    if (BLOCKSIZE > p_nrows) {
      blocksize = p_nrows*p_ncols;
    } else {
      blocksize = BLOCKSIZE*p_ncols;
    }
    double *val_buf = (double*)malloc(blocksize*sizeof(double));
    int *mask_buf = (int*)malloc(blocksize*sizeof(int));
    int jlo = 0;
    int jhi = p_ncols-1;

    // write output in blocks
    int ilo = itask*BLOCKSIZE;
    int ihi = (itask+1)*BLOCKSIZE-1;
    if (ihi >= p_nrows) ihi = p_nrows-1;
    if (ilo<=ihi) {
      lo[0] = ilo; 
      hi[0] = ihi; 
      lo[1] = jlo; 
      hi[1] = jhi; 
      ld = p_ncols;
      NGA_Get(p_data,lo,hi,val_buf,&ld);
      NGA_Get(p_mask,lo,hi,mask_buf,&ld);
      int nrows = ihi-ilo+1;
      int idx;
      vavg.clear();
      vavg2.clear();
      vdiff2.clear();
      for (i=0; i<nrows; i++) {
        idx = i*p_ncols;
        double base = val_buf[idx];
        double avg = 0.0;
        double avg2 = 0.0;
        double diff;
        double diff2 = 0.0;
        int ncnt = 0;
        for (j=0; j<p_ncols; j++) {
          idx = i*p_ncols+j;
          if (mask_buf[idx] >= mval) {
            ncnt++;
            avg += val_buf[idx];
            avg2 += val_buf[idx]*val_buf[idx];
            if (j>0) {
              diff = val_buf[idx] - base;
              diff2 += diff*diff;
            }
          }
        }
        if (ncnt > 0) {
          avg /= ((double)ncnt);
        } else {
          avg = 0.0;
          avg2 = 0.0;
          diff2 = 0.0;
        }
        if (ncnt > 1) {
          avg2 = (avg2-((double)ncnt)*avg*avg)/((double)(ncnt-1));
          diff2 /= ((double)(ncnt-1));
        } else {
          avg2 = 0.0;
          diff2 = 0.0;
        }
        if (avg2 > 0.0) {
          avg2 = sqrt(avg2);
        } else {
          avg2 = 0.0;
        }
        if (diff2 > 0.0) {
          diff2 = sqrt(diff2);
        } else {
          diff2 = 0.0;
        }
        vavg.push_back(avg);
        vavg2.push_back(avg2);
        vdiff2.push_back(diff2);

      }
      // push all results to g_buf
      lo[0] = ilo;
      hi[0] = ihi;
      lo[1] = 0;
      hi[1] = 0;
      ld = 1;
      NGA_Put(g_buf,lo,hi,&vavg[0],&ld);
      lo[1] = 1;
      hi[1] = 1;
      NGA_Put(g_buf,lo,hi,&vavg2[0],&ld);
      lo[1] = 2;
      hi[1] = 2;
      NGA_Put(g_buf,lo,hi,&vdiff2[0],&ld);
    }
    free(val_buf);
    free(mask_buf);
    itask = NGA_Read_inc(g_cnt,&zero,(long)one);
  }
  GA_Pgroup_sync(p_GAgrp);
  // Get data from g_buf and write it to external file
  if (p_me == 0) {
    int ilo = 0;
    int ihi = p_nrows-1;
    char sbuf[128];
    index_set *idx_buf = (index_set*)malloc(p_nrows*sizeof(index_set));
    NGA_Get(p_tags,&ilo,&ihi,idx_buf,&one);
    lo[0] = ilo;
    hi[0] = ihi;
    lo[1] = 0;
    hi[1] = 0;
    ld = 1;
    vavg.resize(p_nrows);
    vavg2.resize(p_nrows);
    vdiff2.resize(p_nrows);
    NGA_Get(g_buf,lo,hi,&vavg[0],&one);
    lo[1] = 1;
    hi[1] = 1;
    NGA_Get(g_buf,lo,hi,&vavg2[0],&one);
    lo[1] = 2;
    hi[1] = 2;
    NGA_Get(g_buf,lo,hi,&vdiff2[0],&one);
    std::ofstream fout;
    fout.open(filename.c_str());
    for (i=0; i<p_nrows; i++) {
      if (flag) {
        if (p_branch_flag) {
          sprintf(sbuf,"%8d %8d %8d %s %16.8e %16.8e %16.8e",idx_buf[i].gidx,
              idx_buf[i].idx1, idx_buf[i].idx2, idx_buf[i].tag, vavg[i], vavg2[i],
              vdiff2[i]);
        } else {
          sprintf(sbuf,"%8d %8d %s %16.8e %16.8e %16.8e",idx_buf[i].gidx,
              idx_buf[i].idx1, idx_buf[i].tag, vavg[i], vavg2[i], vdiff2[i]);
        }
      } else {
        if (p_branch_flag) {
          sprintf(sbuf,"%8d %8d %8d %16.8e %16.8e %16.8e",idx_buf[i].gidx,
              idx_buf[i].idx1, idx_buf[i].idx2, vavg[i], vavg2[i], vdiff2[i]);
        } else {
          sprintf(sbuf,"%8d %8d %16.8e %16.8e %16.8e",idx_buf[i].gidx,
              idx_buf[i].idx1, vavg[i], vavg2[i], vdiff2[i]);
        }
      }
      fout << sbuf << std::endl;
    }
    fout.close();
    free(idx_buf);
  }
  GA_Destroy(g_cnt);
  GA_Destroy(g_buf);
  GA_Pgroup_sync(p_GAgrp);
}

/**
 * Write out file containing Min an Max values in table for each row
 * @param filename name of file containing results
 * @param mval only include values with this mask value or greater
 * @param flag if false, do not include tag ids in output
 */
void stb::writeMinAndMax(std::string filename, int mval, bool flag)
{
  GA_Pgroup_sync(p_GAgrp);
  int zero = 0;
  int one = 1;
  int two = 2;
  int g_cnt = GA_Create_handle();
  NGA_Set_data(g_cnt,one,&one,C_INT);
  NGA_Set_pgroup(g_cnt,p_GAgrp);
  NGA_Allocate(g_cnt);
  GA_Zero(g_cnt);
  int nblock = p_nrows/BLOCKSIZE;
  while (p_nrows > nblock*BLOCKSIZE) {
    nblock++;
  }
  int itask = NGA_Read_inc(g_cnt,&zero,(long)one);
  int dims[2];
  int g_buf = GA_Create_handle();
  dims[0] = p_nrows;
  dims[1] = 5;
  NGA_Set_data(g_buf,two,dims,C_DBL);
  NGA_Set_pgroup(g_buf,p_GAgrp);
  NGA_Allocate(g_buf);
  std::vector<double> vbase, vmin, vmax;
  std::vector<double> idxmin, idxmax;
  int lo[2];
  int hi[2];
  int ld;
  int i, j; 
  while (itask < nblock) {
    // Buffers to hold blocks of data
    int blocksize;
    if (BLOCKSIZE > p_nrows) {
      blocksize = p_nrows*p_ncols;
    } else {
      blocksize = BLOCKSIZE*p_ncols;
    }
    double *val_buf = (double*)malloc(blocksize*sizeof(double));
    int *mask_buf = (int*)malloc(blocksize*sizeof(int));
    int jlo = 0;
    int jhi = p_ncols-1;

    // write output in blocks
    int ilo = itask*BLOCKSIZE;
    int ihi = (itask+1)*BLOCKSIZE-1;
    if (ihi >= p_nrows) ihi = p_nrows-1;
    if (ilo<=ihi) {
      lo[0] = ilo; 
      hi[0] = ihi; 
      lo[1] = jlo; 
      hi[1] = jhi; 
      ld = p_ncols;
      NGA_Get(p_data,lo,hi,val_buf,&ld);
      NGA_Get(p_mask,lo,hi,mask_buf,&ld);
      int nrows = ihi-ilo+1;
      int idx;
      vbase.clear();
      vmin.clear();
      vmax.clear();
      idxmin.clear();
      idxmax.clear();
      for (i=0; i<nrows; i++) {
        idx = i*p_ncols;
        double base = val_buf[idx];
        double min = val_buf[idx];
        double max = val_buf[idx];
        int jmin = 0;
        int jmax = 0;
        for (j=0; j<p_ncols; j++) {
          idx = i*p_ncols+j;
          if (mask_buf[idx] >= mval) {
            if (val_buf[idx] > max) {
              max = val_buf[idx];
              jmax = j;
            }
            if (val_buf[idx] < min) {
              min = val_buf[idx];
              jmin = j;
            }
          }
        }
        vbase.push_back(base);
        vmin.push_back(min);
        vmax.push_back(max);
        idxmin.push_back(static_cast<double>(jmin));
        idxmax.push_back(static_cast<double>(jmax));
      }
      // push all results to g_buf
      lo[0] = ilo;
      hi[0] = ihi;
      lo[1] = 0;
      hi[1] = 0;
      ld = 1;
      NGA_Put(g_buf,lo,hi,&vbase[0],&ld);
      lo[1] = 1;
      hi[1] = 1;
      NGA_Put(g_buf,lo,hi,&vmin[0],&ld);
      lo[1] = 2;
      hi[1] = 2;
      NGA_Put(g_buf,lo,hi,&vmax[0],&ld);
      lo[1] = 3;
      hi[1] = 3;
      NGA_Put(g_buf,lo,hi,&idxmin[0],&ld);
      lo[1] = 4;
      hi[1] = 4;
      NGA_Put(g_buf,lo,hi,&idxmax[0],&ld);
    }
    free(val_buf);
    free(mask_buf);
    itask = NGA_Read_inc(g_cnt,&zero,static_cast<long>(one));
  }
  GA_Pgroup_sync(p_GAgrp);
  // Get data from g_buf and write it to external file
  if (p_me == 0) {
    int ilo = 0;
    int ihi = p_nrows-1;
    char sbuf[256];
    index_set *idx_buf = (index_set*)malloc(p_nrows*sizeof(index_set));
    double *minmax = (double*)malloc(2*p_nrows*sizeof(double));
    NGA_Get(p_tags,&ilo,&ihi,idx_buf,&one);
    lo[0] = ilo;
    hi[0] = ihi;
    lo[1] = 0;
    hi[1] = 0;
    ld = 1;
    vbase.resize(p_nrows);
    vmin.resize(p_nrows);
    vmax.resize(p_nrows);
    idxmin.resize(p_nrows);
    idxmax.resize(p_nrows);
    NGA_Get(g_buf,lo,hi,&vbase[0],&one);
    lo[1] = 1;
    hi[1] = 1;
    NGA_Get(g_buf,lo,hi,&vmin[0],&one);
    lo[1] = 2;
    hi[1] = 2;
    NGA_Get(g_buf,lo,hi,&vmax[0],&one);
    lo[1] = 3;
    hi[1] = 3;
    NGA_Get(g_buf,lo,hi,&idxmin[0],&one);
    lo[1] = 4;
    hi[1] = 4;
    NGA_Get(g_buf,lo,hi,&idxmax[0],&one);
    lo[1] = 0;
    hi[1] = 1;
    NGA_Get(p_bounds,lo,hi,minmax,&two);
    std::ofstream fout;
    fout.open(filename.c_str());
    int idx;
    for (i=0; i<p_nrows; i++) {
      if (flag) {
        if (p_branch_flag) {
          sprintf(sbuf,"%8d %8d %8d %s %16.8e %16.8e %16.8e %16.8e %16.8e",
              idx_buf[i].gidx, idx_buf[i].idx1, idx_buf[i].idx2,
              idx_buf[i].tag, vbase[i], vmin[i], vmax[i],
              vmin[i]-vbase[i], vmax[i]-vbase[i]);
        } else {
          sprintf(sbuf,"%8d %8d %s %16.8e %16.8e %16.8e %16.8e %16.8e",
              idx_buf[i].gidx, idx_buf[i].idx1, idx_buf[i].tag,
              vbase[i], vmin[i], vmax[i], vmin[i]-vbase[i], vmax[i]-vbase[i]);
        }
      } else {
        if (p_branch_flag) {
          sprintf(sbuf,"%8d %8d %8d %16.8e %16.8e %16.8e %16.8e %16.8e",
              idx_buf[i].gidx, idx_buf[i].idx1, idx_buf[i].idx2,
              vbase[i], vmin[i], vmax[i], vmin[i]-vbase[i], vmax[i]-vbase[i]);
        } else {
          sprintf(sbuf,"%8d %8d %16.8e %16.8e %16.8e %16.8e %16.8e",
              idx_buf[i].gidx, idx_buf[i].idx1,
              vbase[i], vmin[i], vmax[i], vmin[i]-vbase[i], vmax[i]-vbase[i]);
        }
      }
      int len = strlen(sbuf);
      char *ptr = sbuf+len;
      if (p_min_bound) {
        idx = i*2;
        sprintf(ptr," %16.8e",minmax[idx]);
      }
      len = strlen(sbuf);
      ptr = sbuf+len;
      if (p_max_bound) {
        idx = i*2+1;
        sprintf(ptr," %16.8e",minmax[idx]);
      }
      len = strlen(sbuf);
      ptr = sbuf+len;
      sprintf(ptr," %8d %8d",static_cast<int>(idxmin[i]),
          static_cast<int>(idxmax[i]));
      fout << sbuf << std::endl;
    }
    fout.close();
    free(minmax);
    free(idx_buf);
  }
  GA_Destroy(g_cnt);
  GA_Destroy(g_buf);
  GA_Pgroup_sync(p_GAgrp);
}

/**
 * Write out file containing number of mask entries at each row that
 * correspond to a given value
 * @param filename name of file containing results
 * @param mval count number of times this mask value occurs
 * @param flag if false, do not include tag ids in output
 */
void stb::writeMaskValueCount(std::string filename, int mval, bool flag)
{
  GA_Pgroup_sync(p_GAgrp);
  int zero = 0;
  int one = 1;
  int two = 2;
  int g_cnt = GA_Create_handle();
  NGA_Set_data(g_cnt,one,&one,C_INT);
  NGA_Set_pgroup(g_cnt,p_GAgrp);
  NGA_Allocate(g_cnt);
  GA_Zero(g_cnt);
  int nblock = p_nrows/BLOCKSIZE;
  while (p_nrows > nblock*BLOCKSIZE) {
    nblock++;
  }
  int itask = NGA_Read_inc(g_cnt,&zero,(long)one);
  int g_buf = GA_Create_handle();
  NGA_Set_data(g_buf,one,&p_nrows,C_INT);
  NGA_Set_pgroup(g_buf,p_GAgrp);
  NGA_Allocate(g_buf);
  std::vector<int> vcnt;
  int lo[2];
  int hi[2];
  int ld;
  int i, j; 
  while (itask < nblock) {
    // Buffers to hold blocks of data
    int blocksize;
    if (BLOCKSIZE > p_nrows) {
      blocksize = p_nrows*p_ncols;
    } else {
      blocksize = BLOCKSIZE*p_ncols;
    }
    int *mask_buf = (int*)malloc(blocksize*sizeof(int));
    int jlo = 0;
    int jhi = p_ncols-1;

    // write output in blocks
    int ilo = itask*BLOCKSIZE;
    int ihi = (itask+1)*BLOCKSIZE-1;
    if (ihi >= p_nrows) ihi = p_nrows-1;
    if (ilo<=ihi) {
      lo[0] = ilo; 
      hi[0] = ihi; 
      lo[1] = jlo; 
      hi[1] = jhi; 
      ld = p_ncols;
      NGA_Get(p_mask,lo,hi,mask_buf,&ld);
      int nrows = ihi-ilo+1;
      int idx;
      vcnt.clear();
      for (i=0; i<nrows; i++) {
        idx = i*p_ncols;
        int icnt = 0;
        for (j=0; j<p_ncols; j++) {
          idx = i*p_ncols+j;
          if (mask_buf[idx] == mval) icnt++;
        }
        vcnt.push_back(icnt);
      }
      // push all results to g_buf
      lo[0] = ilo;
      hi[0] = ihi;
      ld = 1;
      NGA_Put(g_buf,lo,hi,&vcnt[0],&ld);
    }
    free(mask_buf);
    itask = NGA_Read_inc(g_cnt,&zero,(long)one);
  }
  GA_Pgroup_sync(p_GAgrp);
  // Get data from g_buf and write it to external file
  if (p_me == 0) {
    int ilo = 0;
    int ihi = p_nrows-1;
    char sbuf[128];
    index_set *idx_buf = (index_set*)malloc(p_nrows*sizeof(index_set));
    NGA_Get(p_tags,&ilo,&ihi,idx_buf,&one);
    lo[0] = ilo;
    hi[0] = ihi;
    ld = 1;
    vcnt.resize(p_nrows);
    NGA_Get(g_buf,lo,hi,&vcnt[0],&one);
    std::ofstream fout;
    fout.open(filename.c_str());
    for (i=0; i<p_nrows; i++) {
      if (flag) {
        if (p_branch_flag) {
          sprintf(sbuf,"%8d %8d %8d %s %8d", idx_buf[i].gidx,
              idx_buf[i].idx1, idx_buf[i].idx2, idx_buf[i].tag, vcnt[i]);
        } else {
          sprintf(sbuf,"%8d %8d %s %8d", idx_buf[i].gidx,
              idx_buf[i].idx1, idx_buf[i].tag, vcnt[i]);
        }
      } else {
        if (p_branch_flag) {
          sprintf(sbuf,"%8d %8d %8d %8d", idx_buf[i].gidx,
              idx_buf[i].idx1, idx_buf[i].idx2, vcnt[i]);
        } else {
          sprintf(sbuf,"%8d %8d %8d", idx_buf[i].gidx,
              idx_buf[i].idx1, vcnt[i]);
        }
      }
      fout << sbuf << std::endl;
    }
    fout.close();
    free(idx_buf);
  }
  GA_Destroy(g_cnt);
  GA_Destroy(g_buf);
  GA_Pgroup_sync(p_GAgrp);
}

/**
 * Sum up the values in the columns and print the result as a function
 * of column index
 * @param filename name of file containing results
 * @param mval only include values with this mask value or
 greater
 */
void stb::sumColumnValues(std::string filename, int mval)
{
  GA_Pgroup_sync(p_GAgrp);
  int zero = 0;
  int one = 1;
  int two = 2;
  int g_cnt = GA_Create_handle();
  NGA_Set_data(g_cnt,one,&one,C_INT);
  NGA_Set_pgroup(g_cnt, p_GAgrp);
  NGA_Allocate(g_cnt);
  GA_Zero(g_cnt);
  int nblock = p_ncols/BLOCKSIZE;
  while (p_ncols > nblock*BLOCKSIZE) {
    nblock++;
  }
  int itask = NGA_Read_inc(g_cnt,&zero,(long)one);
  int g_buf = GA_Create_handle();
  NGA_Set_data(g_buf,one,&p_ncols,C_DBL);
  NGA_Set_pgroup(g_buf, p_GAgrp);
  NGA_Allocate(g_buf);
  std::vector<double> vsum;
  int lo[2];
  int hi[2];
  int ld;
  int i, j; 
  while (itask < nblock) {
    // Buffers to hold blocks of data
    int blocksize;
    if (BLOCKSIZE > p_ncols) {
      blocksize = p_nrows*p_ncols;
    } else {
      blocksize = BLOCKSIZE*p_nrows;
    }
    double *val_buf = (double*)malloc(blocksize*sizeof(double));
    int *mask_buf = (int*)malloc(blocksize*sizeof(int));
    int ilo = 0;
    int ihi = p_nrows-1;

    // write output in blocks
    int jlo = itask*BLOCKSIZE;
    int jhi = (itask+1)*BLOCKSIZE-1;
    if (jhi >= p_ncols) jhi = p_ncols-1;
    if (jlo<=jhi) {
      lo[0] = ilo; 
      hi[0] = ihi; 
      lo[1] = jlo; 
      hi[1] = jhi; 
      ld = jhi-jlo+1;
      NGA_Get(p_data,lo,hi,val_buf,&ld);
      NGA_Get(p_mask,lo,hi,mask_buf,&ld);
      int ncols = jhi-jlo+1;
      int idx;
      vsum.clear();
      for (j=0; j<ncols; j++) {
        double sum = 0.0;
        for (i=0; i<p_nrows; i++) {
          idx = i*ld+j;
          if (mask_buf[idx] >= mval) {
            sum += val_buf[idx];
          }
        }
        vsum.push_back(sum);
      }
      // push all results to g_buf
      lo[0] = jlo;
      hi[0] = jhi;
      ld = 1;
      NGA_Put(g_buf,lo,hi,&vsum[0],&ld);
    }
    free(mask_buf);
    free(val_buf);
    itask = NGA_Read_inc(g_cnt,&zero,(long)one);
  }
  GA_Pgroup_sync(p_GAgrp);
  // Get data from g_buf and write it to external file
  if (p_me == 0) {
    int ilo = 0;
    int ihi = p_ncols-1;
    char sbuf[128];
    ld = 1;
    vsum.resize(p_ncols);
    NGA_Get(g_buf,&ilo,&ihi,&vsum[0],&one);
    std::ofstream fout;
    fout.open(filename.c_str());
    for (i=0; i<p_ncols; i++) {
      double sum_avg=0.0;
      if (p_nrows > 0) sum_avg = vsum[i]/(static_cast<double>(p_nrows));
      sprintf(sbuf,"%8d %16.8e %16.8e",i+ilo,vsum[i],sum_avg);
      fout << sbuf << std::endl;
    }
    fout.close();
  }
  GA_Destroy(g_cnt);
  GA_Destroy(g_buf);
  GA_Pgroup_sync(p_GAgrp);
}
