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
#include "gridpack/parallel/index_hash.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/environment/environment.hpp"

#define XDIM 100
#define YDIM 100

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
  gridpack::parallel::Communicator world;
  me = world.rank();
  int nprocs;
  nprocs = world.size();

  // Get timer
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();

  // Factor processors into a processor grid
  int ipx, ipy, pdx, pdy;
  factor_grid(nprocs, XDIM, YDIM, &pdx, &pdy);
  if (me == 0) {
    printf("\nProcessor configuration is %d X %d\n",pdx,pdy);
  }
  ipx = me%pdx;
  ipy = (me-ipx)/pdx;

  int ixmin, ixmax, iymin, iymax; // bounds of locally owned nodes
  ixmin = static_cast<int>((static_cast<double>(ipx*XDIM))/(static_cast<double>(pdx)));
  ixmax = static_cast<int>((static_cast<double>((ipx+1)*XDIM))/(static_cast<double>(pdx)))-1;
  iymin = static_cast<int>((static_cast<double>(ipy*YDIM))/(static_cast<double>(pdy)));
  iymax = static_cast<int>((static_cast<double>((ipy+1)*YDIM))/(static_cast<double>(pdy)))-1;

  // Add grid node values to hash table
  int n, ix, iy, nx, ny, i, j;
  nx = ixmax - ixmin + 1;
  ny = iymax - iymin + 1;
  std::vector<std::pair<int,int> > node_pairs;
  for (j=0; j<ny; j++) {
    iy = j + iymin;
    for (i=0; i<nx; i++) {
      ix = i + ixmin;
      n = iy*XDIM + ix;
      node_pairs.push_back(std::pair<int,int>(n,2*n));
    }
  }
  int t_init = timer->createCategory("Initialize Index Hash Object");
  timer->start(t_init);
  gridpack::hash_map::GlobalIndexHashMap hashMap(world);
  timer->stop(t_init);
  int t_keys = timer->createCategory("Add Single Keys");
  timer->start(t_keys);
  hashMap.addPairs(node_pairs);
  timer->stop(t_keys);

  // Add edges to hash table. Start with edges connecting nodes in the
  // i-direction
  int n1, n2;
  std::vector<std::pair<std::pair<int,int>,int> > edge_pairs;
  if (ixmax != XDIM-1) {
    nx = ixmax - ixmin + 1;
  } else {
    nx = ixmax - ixmin;
  }
  ny = iymax - iymin + 1;
  for (j=0; j<ny; j++) {
    iy = j + iymin;
    for (i=0; i<nx; i++) {
      ix = i + ixmin;
      n1 = iy*XDIM+ix;
      n2 = iy*XDIM+ix+1;
      n = iy*(XDIM-1) + ix;
      edge_pairs.push_back(std::pair<std::pair<int,int>,int>
          (std::pair<int,int>(n1,n2),n));
    }
  }
  // Add edges in the j-direction to hash table
  nx = ixmax - ixmin + 1;
  if (iymax != YDIM) {
    ny = iymax - iymin + 1;
  } else {
    ny = iymax - iymin;
  }
  for (j=0; j<ny; j++) {
    iy = j + iymin;
    for (i=0; i<nx; i++) {
      ix = i + ixmin;
      n1 = iy*XDIM+ix;
      n2 = (iy+1)*XDIM+ix;
      n = iy*XDIM + ix + (XDIM-1)*YDIM;
      edge_pairs.push_back(std::pair<std::pair<int,int>,int>
          (std::pair<int,int>(n1,n2),n));
    }
  }
  int t_2key = timer->createCategory("Add Key Pairs");
  timer->start(t_2key);
  hashMap.addPairs(edge_pairs);
  timer->stop(t_2key);

  // Network is complete. Now set up data structures that are to go to each bus
  // and branch
  int nnode = XDIM*YDIM;
  int nedge = (XDIM-1)*YDIM + XDIM*(YDIM-1);
  int half_edge = (XDIM-1)*YDIM;
  double r_me = static_cast<double>(me);
  double r_nnode = static_cast<double>(nnode);
  double r_nedge = static_cast<double>(nedge);
  double r_nproc = static_cast<double>(nprocs);
  int min_node = static_cast<int>(r_me*r_nnode/r_nproc);
  int max_node = static_cast<int>((r_me+1.0)*r_nnode/r_nproc) - 1;
  int min_edge = static_cast<int>(r_me*r_nedge/r_nproc);
  int max_edge = static_cast<int>((r_me+1.0)*r_nedge/r_nproc)-1;

  // Test node values
  std::vector<int> node_keys;
  std::vector<int> node_values;
  for (i=min_node; i<=max_node; i++) {
    node_keys.push_back(i);
  }
  int t_keyval = timer->createCategory("Get Single Key Values");
  timer->start(t_keyval);
  hashMap.getValues(node_keys,node_values);
  timer->stop(t_keyval);
  bool ok = true;
  int nsize = node_keys.size();
  for (i=0; i<nsize; i++) {
    if (node_values[i] != 2*node_keys[i]) {
      printf("p[%d] mismatched nodes 2*key: %d value: %d\n",me,2*node_keys[i],
          node_values[i]);
      ok = false;
    }
  }
  int oks, okr;
  BOOST_CHECK(ok);
  oks = static_cast<int>(ok);
  ierr = MPI_Allreduce(&oks, &okr, 1, MPI_INT, MPI_PROD, mpi_world);
  ok = static_cast<bool>(okr);
  if (me == 0 && ok) {
    printf("\nKey test succeeded\n");
  } else if (me == 0 && !ok) {
    printf("\nKey test failed\n");
  }

  // Test edge values
  std::vector<std::pair<int,int> > edge_keys;
  std::vector<int> edge_values;
  for (i = min_node; i<=max_node; i++) {
    ix = i%XDIM;
    iy = (i-ix)/XDIM;
    if (ix < XDIM-1) {
      n = iy*XDIM + ix + 1;
      edge_keys.push_back(std::pair<int,int>(i,n));
    }
    if (iy < YDIM-1) {
      n = (iy+1)*XDIM + ix;
      edge_keys.push_back(std::pair<int,int>(i,n));
    }
  }
  int t_2keyval = timer->createCategory("Get Key Pair Values");
  timer->start(t_2keyval);
  hashMap.getValues(edge_keys,edge_values);
  timer->stop(t_2keyval);
  nsize = edge_keys.size();
  int ix1, ix2, iy1, iy2;
  ok = true;
  for (i=0; i<nsize; i++) {
    n1 = edge_keys[i].first;
    n2 = edge_keys[i].second;
    ix1 = n1%XDIM;
    iy1 = (n1-ix1)/XDIM;
    ix2 = n2%XDIM;
    iy2 = (n2-ix2)/XDIM;
    if (ix2 == ix1+1) {
      n = iy1*(XDIM-1) + ix1;
      if (n != edge_values[i]) {
        printf("p[%d] mismatched edge values n: %d value: %d\n",me,n,
            edge_values[i]);
        ok = false;
      }
    } else if (iy2 == iy1+1) {
      n = iy1*XDIM + ix1 + (XDIM-1)*YDIM;
      if (n != edge_values[i]) {
        printf("p[%d] mismatched edge values n: %d value: %d\n",me,n,
            edge_values[i]);
        ok = false;
      }
    } else {
      printf("p[%d] Unknown key pair found < %d, %d> <%d, %d>\n",
          me,ix1,iy1,ix2,iy2);
      ok = false;
    }
  }
  BOOST_CHECK(ok);
  oks = static_cast<int>(ok);
  ierr = MPI_Allreduce(&oks, &okr, 1, MPI_INT, MPI_PROD, mpi_world);
  ok = static_cast<bool>(okr);
  if (me == 0 && ok) {
    printf("\nPair key test succeeded\n\n");
  } else if (me == 0 && !ok) {
    printf("\nPair key test failed\n\n");
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
