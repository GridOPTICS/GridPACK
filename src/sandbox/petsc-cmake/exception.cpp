/**
 * @file   exception.cpp
 * @author William A. Perkins
 * @date   2015-06-04 13:51:34 d3g096
 * 
 * @brief  
 * 
 * 
 */


#include <iostream>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>


namespace mpi = boost::mpi;

#include <petscsys.h>

// You gotta love PETSc consistency 

#if PETSC_VERSION_(3,4,0)
#undef PETSC_VERSION_RELEASE
#define PETSC_VERSION_RELEASE 0
#endif

// With PETSc version 3.5, PETSc::Exception was no longer defined.  It
// was replaced with std::runtime_error. 

#if PETSC_VERSION_LT(3,5,0)
#include <petscsys.hh>
#define PETSC_EXCEPTION_TYPE PETSc::Exception
#else
#define PETSC_EXCEPTION_TYPE std::runtime_error
#endif 

void 
doit(const int& me)
{
  try {
    // SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SYS, "bogus system error");
    CHKERRXX(PETSC_ERR_SYS);
    std::cout << me << ": nothing thrown" << std::endl;
  } catch (const PETSC_EXCEPTION_TYPE& e) {
#if PETSC_VERSION_LT(3,5,0)
    std::cout << me << ": caught petsc exception: " << e.message() << std::endl;
#else
    std::cout << me << ": caught petsc exception: " << e.what() << std::endl;
#endif 
  } catch (const std::exception& e) {
    std::cout << me << ": caught exception: " << e.what() << std::endl;
  } catch (...) {
    std::cout << me << ": caught something else" << std::endl;
  }
}  



int 
main(int argc, char* argv[])
{
  // Instantiating a boost::mpi::environment here calls MPI_Init() so
  // PETSc won't
  mpi::environment env(argc, argv);
  mpi::communicator world;

  int nproc = world.size();
  int me = world.rank();

  PetscErrorCode ierr(0);

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);CHKERRXX(ierr);
  // ierr = PetscPopErrorHandler();CHKERRXX(ierr);
  // ierr = PetscPushErrorHandler(&PetscReturnErrorHandler, NULL);CHKERRXX(ierr);

  std::cout << "I am process " << me << " of " << nproc << "." << std::endl;

  world.barrier();

  doit(me);

  ierr = PetscFinalize();CHKERRXX(ierr);

  world.barrier();

  return 0;
}

