// -------------------------------------------------------------
/**
 * @file   bad_destroy.cpp
 * @author William A. Perkins
 * @date   2013-04-23 09:21:01 d3g096
 * 
 * @brief This is a test to investigate what happens if one tries to
 * destroy a PETSc object after PETSc has been finalized. 
 * 
 * 
 */
// -------------------------------------------------------------
#include <iostream>
#include <complex>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/utility.hpp>

namespace mpi = boost::mpi;

#include <petscvec.h>

// -------------------------------------------------------------
//  class VWrap
// -------------------------------------------------------------
class VWrap : private boost::noncopyable {
protected:
  Vec x;

  void init(const mpi::communicator& comm) 
  {
    PetscErrorCode ierr;
    PetscInt lo, hi;
    ierr = VecCreate(comm,&x); // CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, 5*comm.size()); // CHKERRQ(ierr);
    ierr = VecSetFromOptions(x); // CHKERRQ(ierr);

    ierr = VecGetOwnershipRange(x, &lo, &hi);
    for (PetscInt i = lo; i <= hi; ++i) {
      std::complex<double> v(i, 5*comm.size() - i - 1);
      ierr = VecSetValue(x, i, v, INSERT_VALUES);
    }

    ierr = VecAssemblyBegin(x); // CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x); // CHKERRQ(ierr);
  }

public:

  /// Default constructor.
  VWrap(const mpi::communicator& comm)
  {
    init(comm);
  }

  /// Destructor
  ~VWrap(void)
  {
    // Bad things happen (e.g. race condition on RHEL5) if one tries
    // to destroy a PETSc thing after PETSc is finalized.
    PetscErrorCode ierr;
    PetscBool ok;
    ierr = PetscInitialized(&ok);
    if (ok) {
      ierr = VecDestroy(&x);
    }
  }

  void view(void) const
  {
    PetscErrorCode ierr;
    ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD); // CHKERRQ(ierr);
  }

};



int main(int argc, char* argv[])
{
  // Instantiating a boost::mpi::environment here calls MPI_Init() so
  // PETSc won't
  mpi::environment env(argc, argv);
  mpi::communicator world;

  int nproc = world.size();
  int me = world.rank();

  PetscErrorCode ierr(0);

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);

  std::cout << "I am process " << me << " of " << nproc << "." << std::endl;

  world.barrier();
  
  VWrap *v(new VWrap(world));

  v->view();
  
  PetscBool ok;

  delete v;

  ierr = PetscFinalize();CHKERRQ(ierr);

  world.barrier();

  return 0;
}

