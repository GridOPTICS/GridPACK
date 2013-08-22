#include <vector>

#include "mpi.h"
#include <macdecls.h>
#include "gridpack/utilities/complex.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"

#define XDIM 100
#define YDIM 100

class TestBus
  : public gridpack::component::BaseBusComponent {
  public: 

  TestBus(void) {
  }

  ~TestBus(void) {
  }

  bool matrixDiagSize(int *isize, int *jsize) const {
    if (!getReferenceBus()) {
      *isize = 1;
      *jsize = 1;
      return true;
    } else {
      *isize = 0;
      *jsize = 0;
      return false;
    }
  }

  bool matrixDiagValues(gridpack::ComplexType *values) {
    if (!getReferenceBus()) {
      *values = -4.0;
      return true;
    } else {
      return false;
    }
  }

  bool vectorSize(int *isize) const {
    if (!getReferenceBus()) {
      *isize = 1;
      return true;
    } else {
      *isize = 0;
      return false;
    }
  }

  bool vectorValues(gridpack::ComplexType *values) {
    if (!getReferenceBus()) {
      int idx;
      getMatVecIndex(&idx);
      *values = (double)idx;
      return true;
    } else {
      return false;
    }
  }
  
  void setValues(gridpack::ComplexType *values) {
    p_val = *values;
  }

  double getValue(){
    return real(p_val);
  }

  gridpack::ComplexType p_val;
};

class TestBranch
  : public gridpack::component::BaseBranchComponent {
  public: 

  TestBranch(void) {
  }

  ~TestBranch(void) {
  }

  bool matrixForwardSize(int *isize, int *jsize) const {
    if (checkReferenceBus()) {
      *isize = 1;
      *jsize = 1;
      return true;
    } else {
      *isize = 0;
      *jsize = 0;
      return false;
    }
  }

  bool matrixReverseSize(int *isize, int *jsize) const {
    if (checkReferenceBus()) {
      *isize = 1;
      *jsize = 1;
      return true;
    } else {
      *isize = 0;
      *jsize = 0;
      return false;
    }
  }

  bool matrixForwardValues(gridpack::ComplexType *values) {
    if (checkReferenceBus()) {
      *values = 1.0;
      return true;
    } else {
      return false;
    }
  }

  bool matrixReverseValues(gridpack::ComplexType *values) {
    if (checkReferenceBus()) {
      *values = 1.0;
      return true;
    } else {
      return false;
    }
  }

  bool checkReferenceBus() const {
    bool ret = true;
    TestBus *bus1 = dynamic_cast<TestBus*>(getBus1().get());
    TestBus *bus2 = dynamic_cast<TestBus*>(getBus2().get());
    ret = ret && !bus1->getReferenceBus();
    ret = ret && !bus2->getReferenceBus();
    return ret;
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

  if (me == 0) {
    printf("\nTesting FullMatrixMap\n");
  }
  gridpack::mapper::FullMatrixMap<TestNetwork> mMap(network); 
  boost::shared_ptr<gridpack::math::Matrix> M = mMap.mapToMatrix();
  mMap.mapToMatrix(M);

  // Check to see if matrix has correct values
  int one = 1;
  int chk = 0;
  int nbus = network->numBuses();
  int nbranch = network->numBranches();
  int idx,jdx,isize,jsize;
  gridpack::ComplexType v;
  double rv;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      if (network->getBus(i)->matrixDiagSize(&isize,&jsize)) {
        network->getBus(i)->getMatVecIndex(&idx);
        idx--;
        M->get_element(idx,idx,v);
        rv = real(v);
        if (rv != -4.0) {
          printf("p[%d] Diagonal matrix error i: %d j:%d v: %f\n",me,idx,idx,rv);
          chk = 1;
        }
      }
    }
  }
  // Get min and max row indices
  int rlo, rhi; 
  rhi = 0;
  rlo = XDIM*YDIM;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      network->getBus(i)->getMatVecIndex(&idx);
      if (rhi<idx) rhi = idx;
      if (rlo>idx) rlo = idx;
    }
  }
  for (i=0; i<nbranch; i++) {
    if (network->getBranch(i)->matrixForwardSize(&isize,&jsize)) {
      network->getBranch(i)->getMatVecIndices(&idx,&jdx);
      idx--;
      jdx--;
      if (idx >= rlo-1 && idx <= rhi-1) {
        M->get_element(idx,jdx,v);
        rv = real(v);
        if (rv != 1.0) {
          printf("p[%d] Forward matrix error i: %d j:%d v: %f\n",me,idx,jdx,rv);
          chk = 1;
        }
      }
      if (jdx >= rlo-1 && jdx <= rhi-1) {
        M->get_element(jdx,idx,v);
        rv = real(v);
        if (rv != 1.0) {
          printf("p[%d] Reverse matrix error i: %d j:%d v: %f\n",me,jdx,idx,rv);
          chk = 1;
        }
      }
    }
  }
  GA_Igop(&chk,one,"+");
  if (me == 0) {
    if (chk == 0) {
      printf("\nMatrix elements are ok\n");
    } else {
      printf("\nError found in matrix elements\n");
    }
  }

  if (me == 0) {
    printf("\nTesting BusVectorMap\n");
  }
  gridpack::mapper::BusVectorMap<TestNetwork> vMap(network); 
  boost::shared_ptr<gridpack::math::Vector> V = vMap.mapToVector();
  vMap.mapToVector(V);

  // Check to see if vector has correct values
  chk = 0;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      if (network->getBus(i)->vectorSize(&isize)) {
        network->getBus(i)->getMatVecIndex(&idx);
        idx--;
        V->get_element(idx,v);
        rv = real(v);
        if (rv != (double)(idx+1)) {
          printf("p[%d] vector error i: %d v: %f\n",me,idx,rv);
          chk = 1;
        }
      }
    }
  }
  GA_Igop(&chk,one,"+");
  if (me == 0) {
    if (chk == 0) {
      printf("\nVector elements are ok\n");
    } else {
      printf("\nError found in vector elements\n");
    }
  }

  if (me == 0) {
    printf("\nTesting mapToBus\n");
  }

  // Multiply values in vector by a factor of 2
  int lo, hi;
  V->local_index_range(lo,hi);
  for (i=lo; i<hi; i++) {
    V->get_element(i,v);
    v *= 2.0;
    V->set_element(i,v);
  }
  // Push values back onto buses
  vMap.mapToBus(V);

  // Check to see if buses have correct values
  chk = 0;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      if (network->getBus(i)->vectorSize(&isize)) {
        network->getBus(i)->getMatVecIndex(&idx);
        rv = network->getBus(i)->getValue();
        if (rv != (double)(2*idx)) {
          printf("p[%d] Bus error i: %d v: %f expected: %f\n",me,idx,rv,double(2*idx));
          chk = 1;
        }
      }
    }
  }
  GA_Igop(&chk,one,"+");
  if (me == 0) {
    if (chk == 0) {
      printf("\nBus values are ok\n");
    } else {
      printf("\nError found in bus value\n");
    }
  }
}

main (int argc, char **argv) {

  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);
  int me;
  // Initialize Math libraries
  gridpack::math::Initialize();

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &me);
  int nprocs;
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);
  if (me == 0) {
    printf("Testing Mapper Module\n");
    printf("\nTest Network is %d X %d\n",XDIM,YDIM);
  }

  run(me, nprocs);

  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
}
