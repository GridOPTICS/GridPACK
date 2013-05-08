/**
 * @file   exception.cpp
 * @author William A. Perkins
 * @date   2013-05-08 13:09:15 d3g096
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
#include <petscsys.hh>
using namespace PETSc;

void 
doit(const int& me)
{
  try {
    // SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SYS, "bogus system error");
    CHKERRXX(PETSC_ERR_SYS);
    std::cout << me << ": nothing thrown" << std::endl;
  } catch (const PETSc::Exception& e) {
    std::cout << me << ": caught petsc exception: " << e.message() << std::endl;
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

