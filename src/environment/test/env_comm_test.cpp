
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   env_test.cpp
 * @author Bruce Palmer
 * @date   2022-11-03 11:25:58 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include <cstdio>
#include <mpi.h>
#include "gridpack/include/gridpack.hpp"


// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  int global_size, global_rank, color, color_size, color_rank;

  MPI_Init(&argc, &argv);

  MPI_Comm world = MPI_COMM_WORLD;
  MPI_Comm_size(world, &global_size);
  MPI_Comm_rank(world, &global_rank);

  MPI_Comm csplit;

  color = global_rank % 2;

  MPI_Comm_split(world, color, global_rank, &csplit);

  MPI_Comm_size(csplit, &color_size);
  MPI_Comm_rank(csplit, &color_rank);

  printf("MPI %d: Process %d of %d (global: %d of %d)\n",
         color, color_rank, color_size, global_rank, global_size);

  {
    using namespace gridpack;
    using namespace gridpack::parallel;
      
    Environment env(argc,argv,csplit);
    Communicator gp_world;

    printf("GP  %d: Process %d of %d (global: %d of %d)\n",
           color, gp_world.rank(), gp_world.size(),
           global_rank, global_size);

  }

  MPI_Finalize();

    
  return 0;
}
