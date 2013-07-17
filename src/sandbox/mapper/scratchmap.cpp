#include <ga.h>

main () {      // Code sketch for compressing indices assuming that matrix will
               // be built using contributions from both buses and branches.
               // Similar code can be used to construct vectors and matrices
               // that are just built from branches

  int me = GA_Nodeid();
  int nprocs = GA_Nnodes();

  boost::shared_ptr<gridpack::network::BaseNetwork> network;

  int numBuses = network->numBuses();
  int numBranches = network->numBranches();
  int i;
  int numActiveBuses = 0;
  int numActiveBranches = 0;
  for (i = 0; i<numBuses; i++) {
    if (network->getActiveBus(i)) {
      numActiveBuses++;
    }
  }
  for (i = 0; i<numBranches; i++) {
    if (network->getActiveBranch(i)) {
      numActiveBranches++;
    }
  }

  int i, one;
  one = 1;

  // Allocate global arrays to hold matrix block sizes
  int totalBuses = numActiveBuses;
  int totalBranches = numActiveBranches;
  GA_Igop(&totalBuses,one,"+");
  GA_Igop(&totalBranches,one,"+");

  int g_idx = GA_Create_handle();
  GA_Set_data(g_idx, one, &totalBuses, C_INT);
  if (!GA_Allocate(g_idx)) {
    // TODO: some kind of error
  }
  GA_Zero(g_idx);
  int g_jdx = GA_Create_handle();
  GA_Set_data(g_jdx, one, &totalBuses, C_INT);
  if (!GA_Allocate(g_jdx)) {
    // TODO: some kind of error
  }
  GA_Zero(g_jdx);

  // Create arrays to hold matrix size information for buses
  int *ibus_size = new int[numBuses];
  int *jbus_size = new int[numBuses];
  int **ibus_idx = new int*[numBuses];
  for(i = 0; i<numBuses; i++) {
    ibus_idx[i] = new int;
  }

  // Loop over buses and find out which ones contribute to matrix
  int idx, isize, jsize, icnt;
  bool status;
  icnt = 0;
  for (i = 0; i<numBuses; i++) {
    network->getBus(i)->getMatVecIndex(&idx)
    status = network->getBus(i)->getMatrixSize(&isize, &jsize);
    if (status) {
      ibus_size[icnt] = isize;
      jbus_size[icnt] = jsize;
      *(ibus_idx[icnt]) = idx;
      icnt++;
    }
  }
  if (icnt > 0) NGA_Scatter(g_idx, ibus_size, ibus_idx, icnt);
  if (icnt > 0) NGA_Scatter(g_jdx, jbus_size, ibus_idx, icnt);
  GA_Sync();

  // Clean up arrays
  for(i = 0; i<numBuses; i++) {
    delete ibus_idx[i];
  }
  delete [] ibus_idx;
  delete [] ibus_size;
  delete [] jbus_size;
  
  // Create arrays to hold matrix size information for branches
  ibus_size = new int[numBranches];
  jbus_size = new int[numBranches];
  ibus_idx = new int*[numBranches];
  int **jbus_idx = new int*[numBranches];
  for(i = 0; i<numBranches; i++) {
    ibus_idx[i] = new int;
    jbus_idx[i] = new int;
  }

  // Loop over branches and find out which ones contribute to matrix
  icnt = 0;
  int jdx;
  for (i = 0; i<numBranches; i++) {
    network->getBranch(i)->getMatVecIndices(&idx, &jdx);
    status = network->getBranch(i)->getMatrixForwardSize(&isize, &jsize);
    if (status) {
      ibus_size[icnt] = isize;
      jbus_size[icnt] = jsize;
      *(ibus_idx[icnt]) = idx;
      *(jbus_idx[icnt]) = jdx;
      icnt++;
    }
  }
  if (icnt > 0) NGA_Scatter(g_idx, ibus_size, ibus_idx, icnt);
  if (icnt > 0) NGA_Scatter(g_jdx, jbus_size, jbus_idx, icnt);
  icnt = 0;
  for (i = 0; i<numBranches; i++) {
    network->getBranch(i)->getMatVecIndices(&idx, &jdx);
    status = network->getBranch(i)->getMatrixReverseSize(&isize, &jsize);
    if (status) {
      ibus_size[icnt] = isize;
      jbus_size[icnt] = jsize;
      *(ibus_idx[icnt]) = jdx;
      *(jbus_idx[icnt]) = idx;
      icnt++;
    }
  }
  if (icnt > 0) NGA_Scatter(g_idx, ibus_size, ibus_idx, icnt);
  if (icnt > 0) NGA_Scatter(g_jdx, jbus_size, jbus_idx, icnt);
  GA_Sync();

  // Clean up arrays
  for(i = 0; i<numBranches; i++) {
    delete ibus_idx[i];
    delete jbus_idx[i];
  }
  delete [] ibus_idx;
  delete [] jbus_idx;
  delete [] ibus_size;
  delete [] jbus_size;

  // At this point, g_idx and g_jdx should be completely marked with matrix
  // sizes along both the i and j axes. Need to sum across all processors to
  // find actual array dimensions
  int lo, hi, ld;
  int isum; jsum;
  void *ptr;
  NGA_Distribution(g_idx, me, &lo, &hi);
  NGA_Access(g_idx, &lo, &hi, ptr, &ld);
  isum = 0;
  for (i=0; i<hi-lo+1; i++) {
    isum += ((int*)ptr)[i];
  }
  NGA_Release(g_idx, &lo, &hi);
  
  NGA_Distribution(g_jdx, me, &lo, &hi);
  NGA_Access(g_jdx, &lo, &hi, ptr, &ld);
  jsum = 0;
  for (i=0; i<hi-lo+1; i++) {
    jsum += ((int*)ptr)[i];
  }
  NGA_Release(g_jdx, &lo, &hi);

  // Get total array dimensions
  isize = isum;
  GA_Igop(isize,one,"+");
  jsize = jsum;
  GA_Igop(jsize,one,"+");
 
  // Evaluate total offsets for each processor
  int *itmp  = new int[nprocs];
  int *jtmp  = new int[nprocs];

  for (i = 0; i<nprocs; i++) {
    itmp[i] = 0;
    jtmp[i] = 0;
  }
  itmp[me] = isum;
  jtmp[me] = jsum;

  GA_Igop(itmp, nprocs, "+");
  GA_Igop(jtmp, nprocs, "+");

  ioff = 0;
  joff = 0;
  for (i = 1; i <= me; i++) {
    ioff += itmp[i-1];
    joff += jtmp[i-1];
  }
  delete [] itmp;
  delete [] jtmp;

  // Evaluate individual offsets for each processor along i axis
  NGA_Distribution(g_idx, me, &lo, &hi);
  NGA_Access(g_idx, &lo, &hi, ptr, &ld);
  int *ioff_idx = new int[hi-lo+1];
  ioff_idx[0] = ioff;
  for (i=1; i<hi-lo+1; i++) {
    ioff_idx[i] = ioff_idx[i-1] + ((int*)ptr)[i-1];
  }
  NGA_Release(g_idx, &lo, &hi);

  // Create global array that holds all offsets
  int g_ioff = GA_Create_handle();
  GA_set_data(g_ioff, one, &totalBuses, C_INT);
  if (!GA_Allocate(g_ioff)) {
    // TODO: some kind of error
  }
  GA_Zero(g_ioff);
  NGA_Put(g_ioff, &lo, &hi, ioff_idx, &one);
  delete [] ioff_idx;
  
  // Evaluate individual offsets for each processor along j axis
  NGA_Distribution(g_jdx, me, &lo, &hi);
  NGA_Access(g_jdx, &lo, &hi, ptr, &ld);
  int *ioff_jdx = new int[hi-lo+1];
  joff_idx[0] = joff;
  for (i=1; i<hi-lo+1; i++) {
    joff_idx[i] = joff_idx[i-1] + ((int*)ptr)[i-1];
  }
  NGA_Release(g_jdx, &lo, &hi);

  // Create global array that holds all offsets
  int g_joff = GA_Create_handle();
  GA_set_data(g_joff, one, &totalBuses, C_INT);
  if (!GA_Allocate(g_joff)) {
    // TODO: some kind of error
  }
  GA_Zero(g_joff);
  NGA_Put(g_joff, &lo, &hi, joff_idx, &one);
  delete [] joff_idx;

  // At this point, offsets for every bus location in the network are available.
  // Can begin creating matrix
}
