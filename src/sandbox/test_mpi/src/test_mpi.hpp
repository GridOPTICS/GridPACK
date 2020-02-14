/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#ifndef _test_mpi_h
#define _test_mpi_h

#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include <algorithm>
#include <mpi.h>
#include <stdio.h>

// simple object for testing if MPI works for a single processor object
// without using mpirun or mpiexec

class TestSerialMPI {
public:
  
  /**
   * Simple constructure
   */
  TestSerialMPI()
  {
    int argc = 0;
    char **argv;
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("MPI_COMM_WORLD size is %d\n",size);
  }

  void show(const std::vector<std::string> words)
  {
    std::ostream_iterator<std::string> iout(std::cout, ", ");
    std::copy(words.begin(), words.end(), iout);
    std::cout << std::endl;
  }
  
  ~TestSerialMPI()
  {
    MPI_Finalize();
  }
};
#endif
