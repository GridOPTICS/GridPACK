/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hash_test.cpp
 * @author Bruce Palmer
 * @date   June 25, 2014
 * 
 * @brief  A simple test of the GridPACK global index hash functionality
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
#include "gridpack/parallel/index_hash.hpp"

#define NVALUES 10

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  GA_Initialize();
  // Create an artificial scope so that all objects call their destructors
  // before GA_Terminate is called
  if (1) {
    gridpack::parallel::Communicator world;

    int nprocs = world.size();
    int me = world.rank();
    int ntotal = NVALUES*(nprocs-1);
    int i;

    gridpack::hash_map::GlobalIndexHashMap hashmap(world);
    // Create pairs of values and add them to hashmap
    std::vector<std::pair<int, int> > pairs;
    if (me != 0) {
      for (i = 0; i<NVALUES; i++) {
        int ival = (me-1)*NVALUES + i;
        pairs.push_back(std::pair<int, int>(ival,ntotal-1-ival));
      }
    }
    hashmap.addPairs(pairs);

    // Create list of keys corresponding to values on next processor
    int nghbr = (me+1)%nprocs;
    std::vector<int> keys;
    std::vector<int> values;
    if (nghbr != 0) {
      for (i=0; i<NVALUES; i++) {
        keys.push_back((nghbr-1)*NVALUES+i);
      }
    }
    hashmap.getValues(keys, values);

    // Check values to see if they are correct
    bool chk = 1;
    if (nghbr != 0) {
      for (i=0; i<NVALUES; i++) {
        int key = keys[i];
        int ival = values[i];
        if (ival != ntotal-1-key) chk = 0;
      }
      printf("p[%d] value of chk: %d\n",me,chk);
    }
  }

  GA_Terminate();
  return 0;
}

