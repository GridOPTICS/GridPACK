#include <stdio.h>
#include <dsim.hpp>
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

  DSim *dsim = new DSim();

  // Set the configuration file
  dsim->setconfigurationfile(inputfile);

  // Read the data. File names given in config file.
  dsim->readnetworkdatafromconfig();

  // Set up dynamic simulation
  dsim->setup();

  printf("start solving:\n");

  // Solve
  dsim->solve();

  delete(dsim);
  ierr = gridpack_finalize();

  return 0;
}
