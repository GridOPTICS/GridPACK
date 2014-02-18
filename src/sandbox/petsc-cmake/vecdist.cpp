// -------------------------------------------------------------
// file: vecdist.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created February 18, 2014 by William A. Perkins
// Last Change: 2014-02-18 13:55:02 d3g096
// -------------------------------------------------------------

#include <iostream>
#include <cassert>
#include <boost/mpi.hpp>
#include <petsc.h>

// -------------------------------------------------------------
// Vec2Vec
// -------------------------------------------------------------
PetscErrorCode
Vec2Vec(Vec x1, Vec x2)
{
  PetscErrorCode ierr(0);

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)x1,&comm);

  int lo, hi;
  ierr = VecGetOwnershipRange(x1, &lo, &hi); CHKERRQ(ierr);

  IS is1;
  ierr = ISCreateStride(comm, hi-lo, lo, 1, &is1); CHKERRQ(ierr);

  VecScatter scat;
  ierr = VecScatterCreate(x1, is1, x2, is1, &scat); CHKERRQ(ierr);

  ierr = VecScatterBegin(scat, x1, x2, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(scat, x1, x2, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  return ierr;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  int rank(world.rank()), size(world.size());

  assert(size == 2);

  PetscErrorCode ierr(0);
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);

  const int gsize(10);
  const int lsize[2] = { 3*gsize/10, 7*gsize/10 };

  Vec x1, x2;

  ierr = VecCreate(world, &x1); CHKERRQ(ierr);
  ierr = VecSetSizes(x1, lsize[rank], gsize); CHKERRQ(ierr);
  ierr = VecSetFromOptions(x1); CHKERRQ(ierr);

  world.barrier();

  int lo, hi;
  ierr = VecGetOwnershipRange(x1, &lo, &hi); CHKERRQ(ierr);

  for (int i = lo; i < hi; ++i) {
    ierr = VecSetValue(x1, i, i, INSERT_VALUES); CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(x1); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x1); CHKERRQ(ierr);

  ierr = VecView(x1, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  ierr = VecCreate(PETSC_COMM_WORLD, &x2); CHKERRQ(ierr);
  ierr = VecSetSizes(x2, PETSC_DECIDE, gsize); CHKERRQ(ierr);
  ierr = VecSetFromOptions(x2); CHKERRQ(ierr);
  ierr = VecSet(x2, -99); CHKERRQ(ierr);

  ierr = Vec2Vec(x1, x2); CHKERRQ(ierr);

  // IS is1, is2;
  // ierr = ISCreateStride(world, hi-lo, lo, 1, &is1); CHKERRQ(ierr);
  // ierr = VecGetOwnershipRange(x2, &lo, &hi); CHKERRQ(ierr);
  // ierr = ISCreateStride(world, hi-lo, lo, 1, &is2); CHKERRQ(ierr);

  // VecScatter scat;
  // ierr = VecScatterCreate(x1, is1, x2, is1, &scat); CHKERRQ(ierr);

  // ierr = VecScatterBegin(scat, x1, x2, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  // ierr = VecScatterEnd(scat, x1, x2, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = VecView(x2, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}

