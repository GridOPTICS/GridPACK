/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   store_test.cpp
 * @author Bruce Palmer
 * @date   2016-07-12 09:29:08 d3g096
 * 
 * @brief  A simple test of the GridPACK global store module
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
#include "gridpack/parallel/global_store.hpp"
#include "gridpack/environment/environment.hpp"

#define MAX_VEC  1000

#define VEC_LEN  100

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
      printf("Testing GlobalStore on %d processors\n\n",nproc);
    }
    int lo = me*MAX_VEC/nproc;
    int hi = (me+1)*MAX_VEC/nproc-1;
    int i, j;
    gridpack::parallel::GlobalStore<data_type> bank(world);
    // Store vectors in global store object
    for (i=lo; i<=hi; i++) {
      std::vector<data_type> vec;
      for (j=0; j<VEC_LEN+me; j++) {
        data_type item;
        item.ival = j+me;
        item.dval = static_cast<double>(j+me+1);
        vec.push_back(item);
      }
      //printf("p[%d] vec[%d].size: %d\n",me,i,vec.size());
      bank.addVector(i,vec);
    }
    // Upload vectors to global data object
    bank.upload();

    // Check values
    if (me < nproc-1) {
      lo = (me+1)*MAX_VEC/nproc;
      hi = (me+2)*MAX_VEC/nproc-1;
    } else {
      lo = 0;
      hi = MAX_VEC/nproc-1;
    }
    int ichk = me + 1;
    if (me==nproc-1) ichk = 0;
    int chk = 1;
    for (i=lo; i<=hi; i++) {
      std::vector<data_type> vec;
      bank.getVector(i, vec);
      for (j=0; j<VEC_LEN+ichk; j++) {
        bool ok = false;
        if (vec[j].ival == j+ichk &&
            vec[j].dval == static_cast<double>(j+ichk+1)) ok = true;
        if (!ok && j==0) {
          printf("p[%d] Mistake found at (vec[%d])[%d]. Expected ival: %d dval: %d"
              " Actual ival: %f dval: %f\n",me,i,j,j+ichk,vec[j].ival,
              static_cast<double>(j+ichk+1),vec[j].dval);
          chk = 0;
        }
      }
    }
    world.sync();
    world.sum(&chk,1);
    if (chk == nproc && me == 0) {
      printf("Vectors OK\n");
    } else if (chk < nproc && me == 0) {
      printf("Error found in vectors\n");
    }
  }
  return 0;
}

