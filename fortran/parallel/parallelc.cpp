/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// file: parallelc.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created May 14, 2014 by William A. Perkins
// Last Change: 2014-05-14 15:26:07 d3g096
// -------------------------------------------------------------

#include <mpi.h>
#include <ga++.h>

extern "C" void
gridpack_initialize_parallel(int ma_stack, int ma_heap)
{
  int argc(0);
  
  int ierr = MPI_Init(&argc, NULL);
  GA_Initialize();
  MA_init(C_DBL, ma_stack, ma_heap);
}

extern "C" void
gridpack_finalize_parallel(void)
{
  GA_Terminate();
  int ierr = MPI_Finalize();
}
