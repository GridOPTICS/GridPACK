#include <stdio.h>
#include <emt.hpp>
#include <gridpack/include/gridpack.hpp>

int gridpack_initialize(int* argcp, char ***argvp)
{
  int  ierr;
  // Initialize MPI library
  ierr = MPI_Init(argcp,argvp);

  // Initialize GA
  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  // Initialize math libraries
  gridpack::math::Initialize(argcp,argvp);

  return ierr;
}

int gridpack_finalize()
{
  int ierr;
  // Finalize Math libraries
  gridpack::math::Finalize();

  // Terminate GA
  GA_Terminate();

  // Clean up MPI libraries
  ierr = MPI_Finalize();

  return ierr;
}

  
int main(int argc, char **argv)
{
  int ierr;

  ierr = gridpack_initialize(&argc,&argv);
    char inputfile[256];
  if (argc >= 2 && argv[1] != NULL) {
    sprintf(inputfile,"%s",argv[1]);
  } else {
    sprintf(inputfile,"input.xml",argv[1]);
  }

  Emt *emt = new Emt();

  // Set the configuration file
  emt->setconfigurationfile(inputfile);

  // Solve Power flow
  emt->solvepowerflow();

  // Set up dynamic simulation
  emt->setup();

  printf("start solving:\n");

  // Solve
  emt->solve();

  delete(emt);
  ierr = gridpack_finalize();

  return 0;
}
