/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * base_export.hpp
 *
 *  This is a utility class that is used to implement the rest of the export
 *  functionality
 *
 *  Created on: April 4, 2019
 *      Author: Bruce Palmer
 */

#ifndef BASEEXPORT_HPP_
#define BASEEXPORT_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "gridpack/parallel/communicator.hpp"
#include "gridpack/network/base_network.hpp"

#define MAX_STRING_SIZE 1024

namespace gridpack {
namespace expnet {

  // Data structure used to move around text lines
  struct text_line {
    // Global index associated with text line
    int global_idx;
    // Device index associated with text line
    int device_idx;
    // Output text
    char text[MAX_STRING_SIZE];
  }; 

template<class _network>
class BaseExport
{
  public:

    /**
     * Constructor
     */
    explicit BaseExport(gridpack::parallel::Communicator comm) : p_comm(comm)
    {
    }

    /**
     * Destructor
     */
    virtual ~BaseExport(){}


    /**
     * export text to fstream
     * @param fout stream object to export data
     * @param text_data vector of text strings that should be written out
     *                  consecutively, based on values in text_line data
     *                  data structures.
     */
    void writeDataBlock(std::ofstream &fout, std::vector<text_line> &text_data)
    {
      // Loop through all text_line data structures and find highest value of
      // global index
      int i;
      int nsize = text_data.size();
      int gmax = -1;
      for (i=0; i<nsize; i++) {
        if (text_data[i].global_idx > gmax) gmax = text_data[i].global_idx;
      }
      gmax++;
      p_comm.max(&gmax,1);
      if (gmax == 0) return;
      // Create a global array to evaluate offsets
      int grp = p_comm.getGroup();
      int one = 1;
      int ttype = NGA_Register_type(sizeof(text_line));
      int g_offset = NGA_Create_handle();
      if (gmax == 0) return;
      NGA_Set_data(g_offset,one,&gmax,C_INT);
      NGA_Set_pgroup(g_offset,grp);
      NGA_Allocate(g_offset);
      NGA_Zero(g_offset);
      // Evaluate offsets by accumulating number of text strings to each global
      // index
      std::vector<int*> indices;
      indices.resize(nsize);
      std::vector<int> ones;
      ones.resize(nsize);
      for (i=0; i<nsize; i++) {
        indices[i] = &text_data[i].global_idx;
        ones[i] = one;
      }
      if (nsize > 0) {
        NGA_Scatter_acc(g_offset,&ones[0],&indices[0],nsize,&one);
      }
      GA_Pgroup_sync(grp);
      // We now have all the sizes in the g_offset array for each global index.
      // Count the offsets for each global element locally and then get the
      // global offsets
      int lo,hi;
      int me = p_comm.rank();
      int nproc = p_comm.size();
      NGA_Distribution(g_offset,me,&lo,&hi);
      int nelem = hi-lo+1;
      int ld;
      int *ptr;
      if (nelem > 0) {
        NGA_Access(g_offset,&lo,&hi,&ptr,&ld);
      }
      std::vector<int> offsets;
      offsets.resize(nelem);
      if (nelem > 0) offsets[0] = 0;
      int tmp;
      int total = 0;
      if (nelem > 0) total = ptr[0];
      for (i=1; i<nelem; i++) {
        offsets[i] = offsets[i-1]+ptr[i-1];
        total += ptr[i];
      }
      for (i=0; i<nelem; i++) {
        ptr[i] = offsets[i];
      }
      // Get total number of elements across all processors
      std::vector<int> block_sizes;
      block_sizes.resize(nproc);
      for (i=0; i<nproc; i++) block_sizes[i] = 0;
      block_sizes[me] = total;
      p_comm.sum(&block_sizes[0],nproc);
      int loffset = 0;
      for (i=0; i<me; i++) loffset += block_sizes[i];
      total = 0;
      for (i=0; i<nproc; i++) total += block_sizes[i];
      // finish up evaluation of offset array
      for (i=0; i<nelem; i++) ptr[i] = ptr[i]+loffset;
      if (nelem > 0) {
        NGA_Release(g_offset,&lo,&hi);
      }

      // Create global array to hold all text_data
      int g_txt = NGA_Create_handle();
      NGA_Set_pgroup(g_txt,grp);
      NGA_Set_data(g_txt,one,&total,ttype);
      NGA_Allocate(g_txt);
      // Find location of each text_line in g_txt
      std::vector<int> tidx;
      tidx.resize(nsize);
      NGA_Gather(g_offset,&tidx[0],&indices[0],nsize);
      // We now have the offset for the global index. Need to augment it by the
      // device index to get a unique location
      for (i=0; i<nsize; i++) {
        tidx[i] = tidx[i] + text_data[i].device_idx;
        indices[i] = &tidx[i];
      }
      // Scatter all text strings to g_txt
      NGA_Scatter(g_txt,&text_data[0],&indices[0],nsize);
      NGA_Pgroup_sync(grp);

      // Write out strings from rank 0
      if (me == 0) {
        int iproc;
        for (iproc=0; iproc<nproc; iproc++) {
          NGA_Distribution(g_txt,iproc,&lo,&hi);
          std::vector<text_line> lines;
          int nlines = hi-lo+1;
          lines.resize(nlines);
          if (nlines > 0) NGA_Get(g_txt,&lo,&hi,&lines[0],&one);
          for (i=0; i<nlines; i++) {
            fout << lines[i].text;
          }
        }
      }
      // Clean up data
      NGA_Destroy(g_offset);
      NGA_Destroy(g_txt);
      NGA_Deregister_type(ttype);
    }

  private:
    boost::shared_ptr<_network>      p_network;

    gridpack::parallel::Communicator p_comm;
};

} /* namespace export */
} /* namespace gridpack */

#endif /* BASEEXPORT_HPP_ */
