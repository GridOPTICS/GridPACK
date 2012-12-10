#include <iostream>

#include <petsc.h>

#include "dynSim.h"

static char help[] = "Performs dynamic simulation of a power grid system\n \
 -i: input file \n							\
 -a: angle output file \n \
 -v: voltage output file \n";

#undef __FUNCT__
#define __FUNC__ "main"
int main(int argc, char *argv[])
{
  int nproc, me;
  PetscErrorCode ierr;

  PetscInitialize(&argc, &argv, NULL, help);

  MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
  MPI_Comm_rank(PETSC_COMM_WORLD, &me);

  double t0, t1;
  DynSim *dynSim = NULL;

  t0 = MPI_Wtime();
  dynSim = new DynSim(me);

  dynSim->readInputSizes();
  t1 = MPI_Wtime();

  if (me == 0)
    std::cout << "Initialization time: " << t1 - t0 << std::endl;

  t0 = MPI_Wtime();
  dynSim->allocMainData();
  t1 = MPI_Wtime();

  if (me == 0)
    std::cout << "Alloc main data time: " << t1 - t0 << std::endl;

  t0 = MPI_Wtime();
  dynSim->readInputData();
  t1 = MPI_Wtime();

  if (me == 0)
    std::cout << "Read input data time: " << t1 - t0 << std::endl;

  t0 = MPI_Wtime();
  dynSim->buildAdmittanceMatrix();
  t1 = MPI_Wtime();

  if (me == 0)
    std::cout << "Build admittance matrix time: " << t1 - t0 << std::endl;

  delete dynSim;

  ierr = PetscFinalize();

  return 0;
} // main
