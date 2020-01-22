/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#include <vector>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "mpi.h"
#include <macdecls.h>
#include "gridpack/network/base_network.hpp"
#include "gridpack/parser/hash_distr.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/environment/environment.hpp"

#define XDIM 10
#define YDIM 10

#define NUM_TEST_VALS 10;

struct test_data { int idx;
                   int idx1;
                   int idx2;
};

class TestBus
  : public gridpack::component::BaseBusComponent {
  public:

  TestBus(void) {
  }

  ~TestBus(void) {
  }

  std::vector<test_data> p_data;
  std::vector<int*> p_vec;
};

BOOST_CLASS_EXPORT(TestBus)

class TestBranch
  : public gridpack::component::BaseBranchComponent {
  public:

  TestBranch(void) {
  }

  ~TestBranch(void) {
  }

  std::vector<test_data> p_data;
  std::vector<double*> p_vec;
};

BOOST_CLASS_EXPORT(TestBranch)

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

BOOST_AUTO_TEST_SUITE ( TestHashDistribution )

BOOST_AUTO_TEST_CASE( TestHashFunctions )
{

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
  int n1, n2, lx, ly, bridx(0);
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

  // Network is complete. Now set up data structures that are to go to each bus
  // and branch
  int nbus = XDIM*YDIM;
  int nbranch = (XDIM-1)*YDIM + XDIM*(YDIM-1);
  int half_branch = (XDIM-1)*YDIM;
  double r_me = static_cast<double>(me);
  double r_nbus = static_cast<double>(nbus);
  double r_nbranch = static_cast<double>(nbranch);
  double r_nproc = static_cast<double>(nprocs);
  int min_bus = static_cast<int>(r_me*r_nbus/r_nproc);
  int max_bus = static_cast<int>((r_me+1.0)*r_nbus/r_nproc) - 1;
  int min_branch = static_cast<int>(r_me*r_nbranch/r_nproc);
  int max_branch = static_cast<int>((r_me+1.0)*r_nbranch/r_nproc)-1;

  std::vector<int> bus_keys;
  std::vector<test_data> bus_values;
  // Create bus data
  for (i=min_bus; i<=max_bus; i++) {
    ix = i%XDIM;
    iy = (i-ix)/XDIM;
    n = 2*i;
    test_data data;
    data.idx = n;
    n1 = n%3+1;
    for (j=0; j<n1; j++) {
      bus_keys.push_back(n);
      bus_values.push_back(data);
    }
  }
  //Create branch data
  std::vector<std::pair<int,int> > branch_keys;
  std::vector<test_data> branch_values;
  for (i=min_branch; i<=max_branch; i++) {
    test_data data;
    if (i<half_branch) {
      ix = i%(XDIM-1);
      iy = (i-ix)/(XDIM-1);
      n1 = 2*(iy*XDIM+ix);
      n2 = 2*(iy*XDIM+ix+1);
      data.idx1 = n1;
      data.idx2 = n2;
    } else {
      n = i - half_branch;
      ix = n%XDIM;
      iy = (n-ix)/XDIM;
      n1 = 2*(iy*XDIM+ix);
      n2 = 2*((iy+1)*XDIM+ix);
      data.idx1 = n1;
      data.idx2 = n2;
    }
    n = (n1+n2)%3+1;
    for (j=0; j<n; j++) {
      branch_keys.push_back(std::pair<int,int>(n1,n2));
      branch_values.push_back(data);
    }
  }
  // Create HashDistribution object and disperse data to processors that own the
  // corresponding buses/branches
  int t_init = timer->createCategory("Initialize Hash Distribution Object");
  timer->start(t_init);
  gridpack::hash_distr::HashDistribution<TestNetwork, test_data, test_data>
    distr(network);
  timer->stop(t_init);
  int t_bus = timer->createCategory("Distribute Bus Data");
  timer->start(t_bus);
  distr.distributeBusValues(bus_keys, bus_values);
  timer->stop(t_bus);
  std::vector<int> branch_ids;
  int t_branch = timer->createCategory("Distribute Branch Data");
  timer->start(t_branch);
  distr.distributeBranchValues(branch_keys, branch_ids, branch_values);
  timer->stop(t_branch);

  // Copy data to bus and branch objects
  int nsize = bus_keys.size();
  for (i=0; i<nsize; i++) {
    dynamic_cast<TestBus*>(network->getBus(bus_keys[i]).get())->
      p_data.push_back(bus_values[i]);
  }
  nsize = branch_ids.size();
  for (i=0; i<nsize; i++) {
    dynamic_cast<TestBranch*>(network->getBranch(branch_ids[i]).get())->
      p_data.push_back(branch_values[i]);
  }

  // Check that each bus and branch has the correct data
  nbus = network->numBuses();
  nbranch = network->numBranches();
  std::vector<test_data> testData;
  int bus_id, from_bus, to_bus, ndata;
  bool ok = true;
  int oks, okr;
  for (i=0; i<nbus; i++) {
    testData = dynamic_cast<TestBus*>(network->getBus(i).get())->p_data;
    bus_id = network->getOriginalBusIndex(i);
    ndata = testData.size();
    if (ndata != bus_id%3 + 1) {
      ok = false;
      printf("bus %d failed ndata: %d\n",bus_id,ndata);
    }
    for (j=0; j<ndata; j++) {
      if (bus_id != testData[j].idx) {
        ok = false;
        printf("bus %d failed idx: %d\n",bus_id,testData[j].idx);
      }
    }
  }
  BOOST_CHECK(ok);
  oks = static_cast<int>(ok);
  ierr = MPI_Allreduce(&oks, &okr, 1, MPI_INT, MPI_PROD, mpi_world);
  ok = static_cast<bool>(okr);
  if (me == 0 && ok) {
    printf("\nDistribution of buses succeeded\n");
  } else if (me == 0 && !ok) {
    printf("\nDistribution of buses failed\n");
  }
  int ix1,ix2,iy1,iy2;
  ok = true;
  for (i=0; i<nbranch; i++) {
    testData = dynamic_cast<TestBranch*>(network->getBranch(i).get())->p_data;
    network->getOriginalBranchEndpoints(i,&from_bus,&to_bus);
    ndata = testData.size();
    if (ndata != (from_bus+to_bus)%3 + 1) {
      ok = false;
      ix1 = (from_bus/2)%XDIM;
      iy1 = (from_bus/2-ix1)/XDIM;
      ix2 = (to_bus/2)%XDIM;
      iy2 = (to_bus/2-ix1)/XDIM;
      printf("p[%d] branch < %d, %d> location < %d, %d> < %d, %d> failed ndata: %d\n",
          me,from_bus,to_bus,ix1,iy1,ix2,iy2,ndata);
    }
    for (j=0; j<ndata; j++) {
      if (from_bus != testData[j].idx1 || to_bus != testData[j].idx2) {
        ok = false;
        printf("p[%d] branch < %d, %d> failed idx1: %d idx2: %d\n",me,from_bus,to_bus,
            testData[j].idx1,testData[j].idx2);
      }
    }
  }
  BOOST_CHECK(ok);
  oks = static_cast<int>(ok);
  ierr = MPI_Allreduce(&oks, &okr, 1, MPI_INT, MPI_PROD, mpi_world);
  ok = static_cast<bool>(okr);
  if (me == 0 && ok) {
    printf("\nDistribution of branches succeeded\n");
  } else if (me == 0 && !ok) {
    printf("\nDistribution of branches failed\n");
  }

  // Perform tests where all data originates on processor 0
  // Create asymmetric bus data
  nbus = XDIM*YDIM;
  nbranch = (XDIM-1)*YDIM + XDIM*(YDIM-1);
  bus_keys.clear();
  bus_values.clear();
  if (me == 0) {
    for (i=0; i<nbus; i++) {
      ix = i%XDIM;
      iy = (i-ix)/XDIM;
      n = 2*i;
      test_data data;
      data.idx = n;
      n1 = n%3+1;
      for (j=0; j<n1; j++) {
        bus_keys.push_back(n);
        bus_values.push_back(data);
      }
    }
  }
  //Create branch data
  branch_keys.clear();
  branch_values.clear();
  branch_ids.clear();
  if (me == 0) {
    for (i=0; i<nbranch; i++) {
      test_data data;
      if (i<half_branch) {
        ix = i%(XDIM-1);
        iy = (i-ix)/(XDIM-1);
        n1 = 2*(iy*XDIM+ix);
        n2 = 2*(iy*XDIM+ix+1);
        data.idx1 = n1;
        data.idx2 = n2;
      } else {
        n = i - half_branch;
        ix = n%XDIM;
        iy = (n-ix)/XDIM;
        n1 = 2*(iy*XDIM+ix);
        n2 = 2*((iy+1)*XDIM+ix);
        data.idx1 = n1;
        data.idx2 = n2;
      }
      n = (n1+n2)%3+1;
      for (j=0; j<n; j++) {
        branch_keys.push_back(std::pair<int,int>(n1,n2));
        branch_values.push_back(data);
      }
    }
  }
  // Distribute asymmetric bus and branch data
  int t_abus = timer->createCategory("Distribute Asymmetric Bus Data");
  timer->start(t_abus);
  distr.distributeBusValues(bus_keys, bus_values);
  timer->stop(t_abus);
  int t_abranch = timer->createCategory("Distribute Asymmetric Branch Data");
  timer->start(t_abranch);
  distr.distributeBranchValues(branch_keys, branch_ids, branch_values);
  timer->stop(t_abranch);

  // Copy data to bus and branch objects
  nbus = network->numBuses();
  nbranch = network->numBranches();
  for (i=0; i<nbus; i++) {
    dynamic_cast<TestBus*>(network->getBus(i).get())->
      p_data.clear();
  }
  for (i=0; i<nbranch; i++) {
    dynamic_cast<TestBranch*>(network->getBranch(i).get())->
      p_data.clear();
  }
  nsize = bus_keys.size();
  for (i=0; i<nsize; i++) {
    dynamic_cast<TestBus*>(network->getBus(bus_keys[i]).get())->
      p_data.push_back(bus_values[i]);
  }
  nsize = branch_ids.size();
  for (i=0; i<nsize; i++) {
    dynamic_cast<TestBranch*>(network->getBranch(branch_ids[i]).get())->
      p_data.push_back(branch_values[i]);
  }

  // Check that each bus and branch has the correct data
  ok = true;
  for (i=0; i<nbus; i++) {
    testData = dynamic_cast<TestBus*>(network->getBus(i).get())->p_data;
    bus_id = network->getOriginalBusIndex(i);
    ndata = testData.size();
    if (ndata != bus_id%3 + 1) {
      ok = false;
      printf("p[%d] asymmetric bus %d failed ndata: %d\n",me,bus_id,ndata);
    }
    for (j=0; j<ndata; j++) {
      if (bus_id != testData[j].idx) {
        ok = false;
        printf("p[%d] asymmetric bus %d failed idx: %d\n",me,bus_id,testData[j].idx);
      }
    }
  }
  BOOST_CHECK(ok);
  oks = static_cast<int>(ok);
  ierr = MPI_Allreduce(&oks, &okr, 1, MPI_INT, MPI_PROD, mpi_world);
  ok = static_cast<bool>(okr);
  if (me == 0 && ok) {
    printf("\nAsymmetric distribution of buses succeeded\n");
  } else if (me == 0 && !ok) {
    printf("\nAsymmetric distribution of buses failed\n");
  }
  ok = true;
  for (i=0; i<nbranch; i++) {
    testData = dynamic_cast<TestBranch*>(network->getBranch(i).get())->p_data;
    network->getOriginalBranchEndpoints(i,&from_bus,&to_bus);
    ndata = testData.size();
    if (ndata != (from_bus+to_bus)%3 + 1) {
      ok = false;
      ix1 = (from_bus/2)%XDIM;
      iy1 = (from_bus/2-ix1)/XDIM;
      ix2 = (to_bus/2)%XDIM;
      iy2 = (to_bus/2-ix1)/XDIM;
      printf("p[%d] asymmetric branch < %d, %d> location < %d, %d> < %d, %d> failed ndata: %d\n",
          me,from_bus,to_bus,ix1,iy1,ix2,iy2,ndata);
    }
    for (j=0; j<ndata; j++) {
      if (from_bus != testData[j].idx1 || to_bus != testData[j].idx2) {
        ok = false;
        printf("p[%d] asymmetric branch < %d, %d> failed idx1: %d idx2: %d\n",me,from_bus,to_bus,
            testData[j].idx1,testData[j].idx2);
      }
    }
  }
  BOOST_CHECK(ok);
  oks = static_cast<int>(ok);
  ierr = MPI_Allreduce(&oks, &okr, 1, MPI_INT, MPI_PROD, mpi_world);
  ok = static_cast<bool>(okr);
  if (me == 0 && ok) {
    printf("\nAsymmetric distribution of branches succeeded\n");
  } else if (me == 0 && !ok) {
    printf("\nAsymmetric distribution of branches failed\n");
  }

  // Create second HashDistribution object and test dispersing data vectors to
  // processors that own the corresponding buses/branches
  int nvals = NUM_TEST_VALS;
  int k, l;
  // Create bus data
  bus_keys.clear();
  std::vector<int*> bus_vec;
  for (i=min_bus; i<=max_bus; i++) {
    ix = i%XDIM;
    iy = (i-ix)/XDIM;
    n = 2*i;
    n1 = n%3+1;
    for (j=0; j<n1; j++) {
      bus_keys.push_back(n);
      int *idata = new int[nvals];
      for (k=0; k<nvals; k++) {
        idata[k] = n+j+k;
      }
      bus_vec.push_back(idata);
    }
  }
  //Create branch data
  branch_keys.clear();
  std::vector<double*> branch_vec;
  for (i=min_branch; i<=max_branch; i++) {
    if (i<half_branch) {
      ix = i%(XDIM-1);
      iy = (i-ix)/(XDIM-1);
      n1 = 2*(iy*XDIM+ix);
      n2 = 2*(iy*XDIM+ix+1);
    } else {
      n = i - half_branch;
      ix = n%XDIM;
      iy = (n-ix)/XDIM;
      n1 = 2*(iy*XDIM+ix);
      n2 = 2*((iy+1)*XDIM+ix);
    }
    n = (n1+n2)%3+1;
    for (j=0; j<n; j++) {
      double *ddata = new double[nvals];
      for (k=0; k<nvals; k++) {
        ddata[k] = static_cast<double>(n1+n2+j+k);
      }
      branch_keys.push_back(std::pair<int,int>(n1,n2));
      branch_vec.push_back(ddata);
    }
  }
  gridpack::hash_distr::HashDistribution<TestNetwork, int, double>
    distr_v(network);
  distr_v.distributeBusValues(bus_keys, bus_vec, nvals);
  branch_ids.clear();
  distr_v.distributeBranchValues(branch_keys, branch_ids, branch_vec, nvals);

  // Copy data to bus and branch objects
  nsize = bus_keys.size();
  for (i=0; i<nsize; i++) {
    dynamic_cast<TestBus*>(network->getBus(bus_keys[i]).get())->
      p_vec.push_back(bus_vec[i]);
  }
  nsize = branch_ids.size();
  for (i=0; i<nsize; i++) {
    dynamic_cast<TestBranch*>(network->getBranch(branch_ids[i]).get())->
      p_vec.push_back(branch_vec[i]);
  }

  // Check to see if values are correct
  std::vector<int*> test_int;
  ok = true;
  for (i=0; i<nbus; i++) {
    test_int = dynamic_cast<TestBus*>(network->getBus(i).get())->p_vec;
    bus_id = network->getOriginalBusIndex(i);
    ndata = test_int.size();
    if (ndata != bus_id%3 + 1) {
      ok = false;
      printf("p[%d] vector<int> bus %d failed ndata: %d\n",me,bus_id,ndata);
    }
    int *iptr;
    // Data may have gotten shuffled so sort it by increasing values of
    // first integer
    for (k=0; k<ndata; k++) {
      for (l=k; l<ndata; l++) {
        if ((test_int[k])[0] > (test_int[l])[0]) {
          iptr = test_int[k];
          test_int[k] = test_int[l];
          test_int[l] = iptr;
        }
      }
    }
    for (j=0; j<ndata; j++) {
      iptr = test_int[j];
      for (k=0; k<nvals; k++) {
        if (bus_id+j+k != iptr[k]) {
          ok = false;
          printf("p[%d] vector<int> bus %d failed j: %d k: %d ptr: %d\n",
              me,bus_id,j,k,iptr[k]);
        }
      }
    }
  }
  BOOST_CHECK(ok);
  oks = static_cast<int>(ok);
  ierr = MPI_Allreduce(&oks, &okr, 1, MPI_INT, MPI_PROD, mpi_world);
  ok = static_cast<bool>(okr);
  if (me == 0 && ok) {
    printf("\nDistribution of bus vectors succeeded\n");
  } else if (me == 0 && !ok) {
    printf("\nDistribution of bus vectors failed\n");
  }

  std::vector<double*> test_double;
  ok = true;
  for (i=0; i<nbranch; i++) {
    test_double = dynamic_cast<TestBranch*>(network->getBranch(i).get())->p_vec;
    network->getOriginalBranchEndpoints(i,&from_bus,&to_bus);
    ndata = test_double.size();
    if (ndata != (from_bus+to_bus)%3 + 1) {
      ok = false;
      printf("p[%d] vector<double> branch < %d, %d> failed ndata: %d\n",
          me,from_bus,to_bus,ndata);
    }
    double *rptr;
    for (k=0; k<ndata; k++) {
      for (l=k; l<ndata; l++) {
        if ((test_double[k])[0] > (test_double[l])[0]) {
          rptr = test_double[k];
          test_double[k] = test_double[l];
          test_double[l] = rptr;
        }
      }
    }
    for (j=0; j<ndata; j++) {
      rptr = test_double[j];
      for (k=0; k<nvals; k++) {
        if (static_cast<double>(from_bus+to_bus+j+k) != rptr[k]) {
          ok = false;
          printf("p[%d] vector<double> branch < %d, %d> failed j: %d k: %d ptr: %16.2f\n",
              me,from_bus,to_bus,j,k,rptr[k]);
        }
      }
    }
  }
  BOOST_CHECK(ok);
  oks = static_cast<int>(ok);
  ierr = MPI_Allreduce(&oks, &okr, 1, MPI_INT, MPI_PROD, mpi_world);
  ok = static_cast<bool>(okr);
  if (me == 0 && ok) {
    printf("\nDistribution of branch vectors succeeded\n\n");
  } else if (me == 0 && !ok) {
    printf("\nDistribution of branch vectors failed\n\n");
  }

  timer->dump();
}

BOOST_AUTO_TEST_SUITE_END( )

bool init_function(void)
{
  return true;
}

int main (int argc, char **argv) {

  gridpack::Environment env(argc, argv);
  gridpack::parallel::Communicator world;

  int me = world.rank();
  if (me == 0) {
    printf("Testing Network Module\n");
    printf("\nTest Network is %d X %d\n",XDIM,YDIM);
  }

  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  return result;
}
