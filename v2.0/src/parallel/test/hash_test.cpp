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
    int ntotal = NVALUES*nprocs;
    int i;

    // Test list in which all processors contribute
    gridpack::hash_map::GlobalIndexHashMap hashmap(world);
    // Create pairs of values and add them to hashmap
    std::vector<std::pair<int, int> > pairs;
    for (i = 0; i<NVALUES; i++) {
      int ival = me*NVALUES + i;
      pairs.push_back(std::pair<int, int>(ival,ntotal-1-ival));
    }
    hashmap.addPairs(pairs);

    // Create list of keys corresponding to values on next processor
    int nghbr = (me+1)%nprocs;
    std::vector<int> keys;
    std::vector<int> values;
    for (i=0; i<NVALUES; i++) {
      keys.push_back(nghbr*NVALUES+i);
    }
    hashmap.getValues(keys, values);

    // Check values to see if they are correct
    int schk = 1;
    for (i=0; i<NVALUES; i++) {
      int key = keys[i];
      int ival = values[i];
      if (ival != ntotal-1-key) schk = 0;
    }
    int rchk;
    int ierr = MPI_Allreduce(&schk, &rchk, 1, MPI_INT, MPI_SUM, world);
    if (me == 0 && rchk == nprocs) {
      printf("\nSymmetric list passed\n");
    } else if (me == 0) {
      printf("\nSymmetric list failed\n");
    }

    // Test list in which one processor does not contribute
    ntotal = NVALUES*(nprocs-1);
    // Create pairs of values and add them to hashmap
    pairs.clear();
    if (me != 0) {
      for (i = 0; i<NVALUES; i++) {
        int ival = (me-1)*NVALUES + i;
        pairs.push_back(std::pair<int, int>(ival,ntotal-1-ival));
      }
    }
    hashmap.addPairs(pairs);

    // Create list of keys corresponding to values on next processor
    nghbr = (me+1)%nprocs;
    keys.clear();
    values.clear();
    if (nghbr != 0) {
      for (i=0; i<NVALUES; i++) {
        keys.push_back((nghbr-1)*NVALUES+i);
      }
    }
    hashmap.getValues(keys, values);

    // Check values to see if they are correct
    schk = 1;
    if (nghbr != 0) {
      for (i=0; i<NVALUES; i++) {
        int key = keys[i];
        int ival = values[i];
        if (ival != ntotal-1-key) schk = 0;
      }
    }
    ierr = MPI_Allreduce(&schk, &rchk, 1, MPI_INT, MPI_SUM, world);
    if (me == 0 && rchk == nprocs) {
      printf("\nAsymmetric list passed\n");
    } else if (me == 0) {
      printf("\nAsymmetric list failed\n");
    }

    // Repeat test for pair keys list

    // Test list in which all processors contribute
    // Create pairs of values and add them to hashmap
    std::vector<std::pair<std::pair<int,int>, int> > pairs2;
    for (i = 0; i<NVALUES; i++) {
      int ival = me*NVALUES + i;
      std::pair<int,int> key = std::pair<int,int>(ival,ival+1);
      pairs2.push_back(std::pair<std::pair<int,int>, int>(key,ntotal-1-ival));
    }
    hashmap.addPairs(pairs2);

    // Create list of keys corresponding to values on next processor
    nghbr = (me+1)%nprocs;
    std::vector<std::pair<int,int> > keys2;
    values.clear();
    for (i=0; i<NVALUES; i++) {
      int ikey = nghbr*NVALUES+i;
      std::pair<int,int> key = std::pair<int,int>(ikey,ikey+1);
      keys2.push_back(key);
    }
    hashmap.getValues(keys2, values);

    // Check values to see if they are correct
    schk = 1;
    for (i=0; i<NVALUES; i++) {
      int key = keys2[i].first;
      int ival = values[i];
      if (ival != ntotal-1-key) schk = 0;
    }
    rchk;
    ierr = MPI_Allreduce(&schk, &rchk, 1, MPI_INT, MPI_SUM, world);
    if (me == 0 && rchk == nprocs) {
      printf("\nSymmetric pair list passed\n");
    } else if (me == 0) {
      printf("\nSymmetric pair list failed\n");
    }

    // Test list in which one processor does not contribute
    ntotal = NVALUES*(nprocs-1);
    // Create pairs of values and add them to hashmap
    pairs2.clear();
    if (me != 0) {
      for (i = 0; i<NVALUES; i++) {
        int ival = (me-1)*NVALUES + i;
        std::pair<int,int> key = std::pair<int,int>(ival,ival+1);
        pairs2.push_back(std::pair<std::pair<int,int>, int>(key,ntotal-1-ival));
      }
    }
    hashmap.addPairs(pairs2);

    // Create list of keys corresponding to values on next processor
    nghbr = (me+1)%nprocs;
    keys2.clear();
    values.clear();
    if (nghbr != 0) {
      for (i=0; i<NVALUES; i++) {
        int ikey = (nghbr-1)*NVALUES+i;
        keys2.push_back(std::pair<int,int>(ikey,ikey+1));
      }
    }
    hashmap.getValues(keys2, values);

    // Check values to see if they are correct
    schk = 1;
    if (nghbr != 0) {
      for (i=0; i<NVALUES; i++) {
        int key = keys2[i].first;
        int ival = values[i];
        if (ival != ntotal-1-key) schk = 0;
      }
    }
    ierr = MPI_Allreduce(&schk, &rchk, 1, MPI_INT, MPI_SUM, world);
    if (me == 0 && rchk == nprocs) {
      printf("\nAsymmetric pair list passed\n");
    } else if (me == 0) {
      printf("\nAsymmetric pair list failed\n");
    }
  }

  GA_Terminate();
  return 0;
}

