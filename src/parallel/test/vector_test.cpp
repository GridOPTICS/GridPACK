/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   vector_test.cpp
 * @author Bruce Palmer
 * @date   2017-07-6 09:29:08 d3g096
 * 
 * @brief  A simple test of the GridPACK global vector module
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------

#include <iostream>
#include <ga.h>
#include "gridpack/parallel/parallel.hpp"
#include "gridpack/parallel/global_vector.hpp"
#include "gridpack/environment/environment.hpp"


#define VEC_LEN 1000

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------

typedef struct {int ival;
                double dval;
} data_type;

int
main(int argc, char **argv)
{
  gridpack::Environment env(argc, argv);
  // Create an artificial scope so that all objects call their destructors
  // before GA_Terminate is called
  if (1) {
    gridpack::parallel::Communicator world;
    int me = world.rank();
    int nproc = world.size();
    if (me == 0) {
      printf("Testing GlobalVector on %d processors\n\n",nproc);
    }
    int lo, hi;
    if (me < nproc-1) {
      lo = (me+1)*VEC_LEN;
      hi = (me+2)*VEC_LEN-1;
    } else {
      lo = 0;
      hi = VEC_LEN-1;
    }
    if (me == 0) {
      printf("Total vector length is %d\n\n",nproc*VEC_LEN);
    }
    int i, j;
    gridpack::parallel::GlobalVector<data_type> bank(world);
    // Store vectors in vector store object
    std::vector<data_type> vec;
    std::vector<int> idx;
    for (i=lo; i<=hi; i++) {
      data_type item;
      item.ival = i;
      item.dval = static_cast<double>(i+1);
      vec.push_back(item);
      idx.push_back(i);
    }
    bank.addElements(idx,vec);
    // Upload vectors to global data object
    bank.upload();

    // Check values
    int maxdim = nproc*VEC_LEN;
    vec.clear();
    idx.clear();
    i = me;
    while (i<maxdim) {
      idx.push_back(i);
      i += nproc;
    }
    int chk = 1;
    bank.getData(idx, vec);
    for (i=0; i<idx.size(); i++) {
      if (vec[i].ival != idx[i] ||
          vec[i].dval != static_cast<double>(idx[i]+1)) {
        printf("p[%d] Mismatch found in element %d: ival(e): %d ival(a): %d\n",
            me,idx[i],idx[i],vec[i].ival);
        chk = 0;
      }
    }
    world.sync();
    world.sum(&chk,1);
    if (chk == nproc && me == 0) {
      printf("Global vector OK\n");
    } else if (chk < nproc && me == 0) {
      printf("Error found in global vector\n");
    }
    // test getAllData command
    bank.getAllData(vec);
    if (vec.size() != maxdim) {
      if (me == 0) {
        printf("getAllData command does not return all elements\n");
      }
    }
    chk = 1;
    for (i=0; i<idx.size(); i++) {
      if (vec[i].ival != i ||
          vec[i].dval != static_cast<double>(i+1)) {
        printf("p[%d] Mismatch found in element %d: ival(e): %d ival(a): %d\n",
            me,idx[i],idx[i],vec[i].ival);
        chk = 0;
      }
    }
    world.sync();
    world.sum(&chk,1);
    if (chk == nproc && me == 0) {
      printf("Global vector OK for getAllData\n");
    } else if (chk < nproc && me == 0) {
      printf("Error found in global vector for getAllData\n");
    }
    // check resetElements and reload
    vec.clear();
    idx.clear();
    i = me;
    while (i<maxdim) {
      data_type item;
      item.ival = 2*i;
      item.dval = static_cast<double>(2*i+1);
      vec.push_back(item);
      idx.push_back(i);
      i += nproc;
    }
    bank.resetElements(idx,vec);
    bank.reload();
    bank.getAllData(vec);
    chk = 1;
    for (i=0; i<idx.size(); i++) {
      if (vec[i].ival != 2*i ||
          vec[i].dval != static_cast<double>(2*i+1)) {
        printf("p[%d] Mismatch found in element %d: ival(e): %d ival(a): %d\n",
            me,idx[i],idx[i],vec[i].ival);
        chk = 0;
      }
    }
    world.sync();
    world.sum(&chk,1);
    if (chk == nproc && me == 0) {
      printf("Global vector OK for resetElements and reload\n");
    } else if (chk < nproc && me == 0) {
      printf("Error found in global vector for resetElements and reload\n");
    }
  }
  return 0;
}

