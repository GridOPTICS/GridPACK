#include <vector>

#include "mpi.h"
#include <macdecls.h>
#include "gridpack/utilities/complex.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/serial_io/serial_io.hpp"

#define XDIM 10
#define YDIM 10

class TestBus
  : public gridpack::component::BaseBusComponent {
  public: 

  TestBus(void) {
  }

  ~TestBus(void) {
  }

  bool serialWrite(char *string, const char *signal) {
    sprintf(string,"  Bus: %4d      %4d\n",getOriginalIndex(),
        getGlobalIndex());
    return true;
  }
};

class TestBranch
  : public gridpack::component::BaseBranchComponent {
  public: 

  TestBranch(void) {
  }

  ~TestBranch(void) {
  }
  bool serialWrite(char *string, const char *signal) {
    sprintf(string,"  Branch: %4d      %4d      %4d      %4d\n",
        getBus1OriginalIndex(),
        getBus2OriginalIndex(),
        getBus1GlobalIndex(),
        getBus2GlobalIndex());
    return true;
  }
};

void factor_grid(int nproc, int xsize, int ysize, int *pdx, int *pdy)
{
  int i,j,it,ip,ifac,pmax,prime[1000], chk;
  int idx, idy;
  double xtmp, ytmp;
  int fac[1000];

  ip = nproc;

  /* Factor nproc completely. First, find all prime numbers less than or equal
   * to nproc. */
  pmax = 0;
  for (i=2; i<=nproc; i++) {
    chk = 0;
    for (j=0; j<pmax; j++) {
      if (i%prime[j] == 0) {
        chk = 1;
        break;
      }
    }
    if (!chk) {
      prime[pmax] = i;
      pmax++;
    }
  }

  /* Find all prime factors of nproc */
  ifac = 0;
  for (i=0; i<pmax; i++) {
    while(ip%prime[i] == 0 && ip > 1) {
      ifac++;
      fac[ifac] = prime[i];
      ip = ip / prime[i];
    }
  }
  /* Find three factors of nproc such that the resulting three dimensions of
     the simulation cell on each processor are as close as possible to being
     the same size */
  xtmp = xsize;
  ytmp = ysize;
  idx = 1;
  idy = 1;
  for (i=ifac; i>0; i--) {
    if (xtmp >= ytmp) {
      idx = fac[i]*idx;
      xtmp /= fac[i];
    } else {
      idy = fac[i]*idy;
      ytmp /= fac[i];
    }
  }
  *pdx = idx;
  *pdy = idy;
}

typedef gridpack::network::BaseNetwork<TestBus, TestBranch> TestNetwork;

void run (const int &me, const int &nprocs)
{
  // Create network
  gridpack::parallel::Communicator world;
  boost::shared_ptr<TestNetwork> network(new TestNetwork(world));

  // Factor processors into a processor grid
  int ipx, ipy, pdx, pdy;
  factor_grid(nprocs, XDIM, YDIM, &pdx, &pdy);
  if (me == 0) {
    printf("\nProcessor configuration is %d X %d\n",pdx,pdy);
  }
  ipx = me%pdx;
  ipy = (me-ipx)/pdx;

  int ixmin, ixmax, iymin, iymax; // bounds of locally owned nodes
  int iaxmin, iaxmax, iaymin, iaymax; // bounds of locally held nodes
  ixmin = static_cast<int>((static_cast<double>(ipx*XDIM))/(static_cast<double>(pdx)));
  ixmax = static_cast<int>((static_cast<double>((ipx+1)*XDIM))/(static_cast<double>(pdx)))-1;
  iymin = static_cast<int>((static_cast<double>(ipy*YDIM))/(static_cast<double>(pdy)));
  iymax = static_cast<int>((static_cast<double>((ipy+1)*YDIM))/(static_cast<double>(pdy)))-1;

  iaxmin = ixmin - 1;
  if (ixmin == 0) iaxmin = 0;
  iaxmax = ixmax + 1;
  if (ixmax == XDIM - 1) iaxmax = XDIM - 1;

  iaymin = iymin - 1;
  if (iymin == 0) iaymin = 0;
  iaymax = iymax + 1;
  if (iymax == YDIM - 1) iaymax = YDIM - 1;

  // Add buses to network
  int ncnt, n, ix, iy, nx, ny, i, j;
  ncnt = 0;
  nx = iaxmax - iaxmin + 1;
  ny = iaymax - iaymin + 1;
  for (j=0; j<ny; j++) {
    iy = j + iaymin;
    for (i=0; i<nx; i++) {
      ix = i + iaxmin;
      n = iy*XDIM + ix;
      n = 2*n;  // Provide original index that is not equal to global index 
      network->addBus(n);
      // Set active flag for network buses
      if (ix >= ixmin && ix <= ixmax && iy >= iymin && iy <= iymax) {
        network->setActiveBus(ncnt, true);
      } else {
        network->setActiveBus(ncnt, false);
      }
      n = n/2;
      network->setGlobalBusIndex(ncnt, n);
      if (ix == 0 && iy == 0) {
        network->setReferenceBus(ncnt);
      }
      ncnt++;
    }
  }

  // Add branches to network. Start with branches connecting buses in the
  // i-direction
  int n1, n2, lx, ly;
  ncnt = 0;
  nx = iaxmax - iaxmin;
  ny = iymax - iymin + 1;
  for (j=0; j<ny; j++) {
    iy = j + iymin;
    for (i=0; i<nx; i++) {
      ix = i + iaxmin;
      n1 = iy*XDIM+ix;
      n1 = 2*n1;
      n2 = iy*XDIM+ix+1;
      n2 = 2*n2;
      network->addBranch(n1, n2);
      n1 = n1/2;
      n2 = n2/2;
      network->setGlobalBusIndex1(ncnt, n1);
      network->setGlobalBusIndex2(ncnt, n2);
      n = iy*(XDIM-1) + ix;
      network->setGlobalBranchIndex(ncnt,n);
      // Figure out local indices of buses
      lx = ix - iaxmin;
      ly = iy - iaymin;
      n1 = ly*(iaxmax-iaxmin+1) + lx;
      n2 = ly*(iaxmax-iaxmin+1) + lx + 1;
      network->setLocalBusIndex1(ncnt,n1);
      network->setLocalBusIndex2(ncnt,n2);
      // Determine which branches are locally held. Use the rule that if bus 1
      // is local, then branch belongs to this processor
      if (ix >= ixmin && ix <= ixmax && iy >= iymin && iy <= iymax) {
        network->setActiveBranch(ncnt, true);
      } else {
        network->setActiveBranch(ncnt, false);
      }
      ncnt++;
    }
  }
  // Add branches connecting buses in the j-direction
  nx = ixmax - ixmin + 1;
  ny = iaymax - iaymin;
  for (j=0; j<ny; j++) {
    iy = j + iaymin;
    for (i=0; i<nx; i++) {
      ix = i + ixmin;
      n1 = iy*XDIM+ix;
      n1 = 2*n1;
      n2 = (iy+1)*XDIM+ix;
      n2 = 2*n2;
      network->addBranch(n1, n2);
      n1 = n1/2;
      n2 = n2/2;
      network->setGlobalBusIndex1(ncnt, n1);
      network->setGlobalBusIndex2(ncnt, n2);
      n = iy*XDIM + ix + (XDIM-1)*YDIM;
      network->setGlobalBranchIndex(ncnt,n);
      // Figure out local indices of buses
      lx = ix - iaxmin;
      ly = iy - iaymin;
      n1 = ly*(iaxmax-iaxmin+1) + lx;
      n2 = (ly+1)*(iaxmax-iaxmin+1) + lx;
      network->setLocalBusIndex1(ncnt,n1);
      network->setLocalBusIndex2(ncnt,n2);
      // Determine which branches are locally held. Use the rule that if bus 1
      // is local, then branch belongs to this processor
      if (ix >= ixmin && ix <= ixmax && iy >= iymin && iy <= iymax) {
        network->setActiveBranch(ncnt, true);
      } else {
        network->setActiveBranch(ncnt, false);
      }
      ncnt++;
    }
  }

  // Set up some remaining properties of network and network components so that
  // matrix-vector interface is ready to go
  gridpack::factory::BaseFactory<TestNetwork> factory(network);
  factory.setComponents();

  // Create serial bus IO object
  gridpack::serial_io::SerialBusIO<TestNetwork> busIO(128,network);
  busIO.header("\n  Bus Properties\n");
  busIO.header("\n       Original   Global\n");
  busIO.write();

  // Create serial branch IO object
  gridpack::serial_io::SerialBranchIO<TestNetwork> branchIO(128,network);
  branchIO.header("\n  Branch Properties\n");
  branchIO.header("\n         Original1  Original2  Global1 Global2\n");
  branchIO.write();
}

main (int argc, char **argv) {

  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);
  int me;

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &me);
  int nprocs;
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);
  if (me == 0) {
    printf("Testing Serial IO Module\n");
    printf("\nTest Network is %d X %d\n",XDIM,YDIM);
  }

  run(me, nprocs);

  GA_Terminate();

  // Clean up MPI libraries
  ierr = MPI_Finalize();
}
