/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#include <vector>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#include "mpi.h"
#include <macdecls.h>
#include "gridpack/environment/environment.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/parser/bus_table.hpp"
#include "gridpack/timer/coarse_timer.hpp"

#define XDIM 10
#define YDIM 10

class TestBus
  : public gridpack::component::BaseBusComponent {
  public:

  TestBus(void) {
  }

  ~TestBus(void) {
  }
};

class TestBranch
  : public gridpack::component::BaseBranchComponent {
  public:

  TestBranch(void) {
  }

  ~TestBranch(void) {
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

int main(int argc, char **argv)
{

  gridpack::Environment env(argc,argv);
  gridpack::parallel::Communicator comm;
  MPI_Comm mpi_world = static_cast<MPI_Comm>(comm);
  int ierr;
  int me;
  ierr = MPI_Comm_rank(mpi_world, &me);
  int nprocs;
  ierr = MPI_Comm_size(mpi_world, &nprocs);

  // Get timer
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();

  // Create network
  gridpack::parallel::Communicator world;
  typedef gridpack::network::BaseNetwork<TestBus, TestBranch> TestNetwork;
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
      n = n+1;  // Provide original index that is not equal to global index 
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
  int n1, n2, lx, ly, bridx(0);
  ncnt = 0;
  nx = iaxmax - iaxmin;
  ny = iymax - iymin + 1;
  for (j=0; j<ny; j++) {
    iy = j + iymin;
    for (i=0; i<nx; i++) {
      ix = i + iaxmin;
      n1 = iy*XDIM+ix;
      n1 = n1+1;
      n2 = iy*XDIM+ix+1;
      n2 = n2+1;
      network->addBranch(n1, n2);
      n1 = n1-1;
      n2 = n2-1;
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
      n1 = n1+1;
      n2 = (iy+1)*XDIM+ix;
      n2 = n2+1;
      network->addBranch(n1, n2);
      n1 = n1-1;
      n2 = n2-1;
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

  // Network is complete. Now read in table data.
  gridpack::bus_table::BusTable<TestNetwork> table(network);
  std::string file = "table.dat";
  table.readTable(file);
  // Check headers
  std::vector<std::string> headers;
  table.getHeaders(headers);
  if (headers.size() != 2 && me == 0) {
    printf("\nNumber of headers found is incorrect: %d\n",static_cast<int>(headers.size()));
  } else if (me == 0) {
    printf("\nNumber of headers is correct\n");
  }
  if (headers[0] != "1_2_3_4_5_6_7_8_9" && me == 0) {
    printf("\nFirst header is incorrect: (%s)\n",headers[0].c_str());
  } else if (me == 0) {
    printf("\nFirst header is correct\n");
  }
  if (headers[1] != "1 2 3 4 5 6 7 8 9" && me == 0) {
    printf("\nSecond header is incorrect: (%s)\n",headers[1].c_str());
  } else if (me == 0) {
    printf("\nSecond header is correct\n");
  }
  std::vector<int> indices;
  std::vector<std::string> tags;
  std::vector<double> values;
  table.getLocalIndices(indices);
  table.getTags(tags);
  bool ok = true;
  for (i=0; i<100; i++) {
    table.getValues(i,values);
    // check to see if values are correct
    for (j=0; j<values.size(); j++) {
      n = indices[j];
      n1 = network->getOriginalBusIndex(n);
      if (values[j] != static_cast<double>((n1-1)*100+i)) {
        printf("Error for bus %d\n",n1);
        ok = false;
      }
    }
  }
  if (ok && me == 0) {
    printf("\nTable values are ok\n");
  } else if (!ok) {
    printf("\nTable values incorrect\n");
  }
}
