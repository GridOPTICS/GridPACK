/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#include "mpi.h"
#include <vector>
#include <macdecls.h>
#include "gridpack/utilities/complex.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/gen_matrix_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/gen_vector_map.hpp"
#include "gridpack/mapper/gen_slab_map.hpp"

#define XDIM 100
#define YDIM 100
#define NSLAB 20

class TestBus
  : public gridpack::component::BaseBusComponent {
  public: 

  TestBus(void) {
    p_vals = new double*[2];
    p_vals[0] = new double[NSLAB];
    p_vals[1] = new double[NSLAB];
  }

  ~TestBus(void) {
    delete [] p_vals[0];
    delete [] p_vals[1];
    delete [] p_vals;
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

  bool vectorValues(gridpack::RealType *values) {
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

  void setValues(gridpack::RealType *values) {
    p_rval = *values;
  }

  double getValue() {
    return real(p_val);
  }

  double getRealValue() {
    return p_rval;
  }

  int matrixNumRows() const {
    return 1;
  }

  int matrixNumCols() const {
    return 1;
  }

  void matrixSetRowIndex(int irow, int idx) {
    int midx;
    getMatVecIndex(&midx);
    p_row_idx = idx;
  }

  void matrixSetColIndex(int icol, int idx) {
    p_col_idx = idx;
  }

  int matrixGetRowIndex(int idx) {
    return p_row_idx;
  }

  int matrixGetColIndex(int idx) {
    return p_col_idx;
  }

  int matrixNumValues() const {
    std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > nghbrs;
    getNeighborBranches(nghbrs);
    return nghbrs.size()+1;
  }

  void matrixGetValues(gridpack::ComplexType *values, int *rows, int *cols) {
    int i, im, jm, nsize, idx, jdx1, jdx2;
    std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > nghbrs;
    getNeighborBranches(nghbrs);
    nsize = nghbrs.size();
    idx = getGlobalIndex();
    im = matrixGetRowIndex(0);
    for (i=0; i<nsize; i++) {
      jm = nghbrs[i]->matrixGetColIndex(0);
      rows[i] = im;
      cols[i] = jm;
      jdx1 = dynamic_cast<gridpack::component::BaseBranchComponent*>(nghbrs[i].get())
        ->getBus1GlobalIndex();
      jdx2 = dynamic_cast<gridpack::component::BaseBranchComponent*>(nghbrs[i].get())
        ->getBus2GlobalIndex();
      if (idx == jdx1) {
        values[i] = static_cast<gridpack::ComplexType>(idx + jdx2);
      } else {
        values[i] = static_cast<gridpack::ComplexType>(idx + jdx1);
      }
    }
    values[nsize] = idx;
    rows[nsize] = im;
    jm = matrixGetColIndex(0);
    cols[nsize] = jm;
  }

  int vectorNumElements() const
  {
     return 2;
  }

  void vectorSetElementIndex(int ielem, int idx)
  {
    if (ielem == 0) {
      p_vec_idx1 = idx;
    } else {
      p_vec_idx2 = idx;
    }
  }

  void vectorGetElementIndices(int *idx)
  {
    idx[0] = p_vec_idx1;
    idx[1] = p_vec_idx2;
  }

  void vectorGetElementValues(gridpack::ComplexType *values, int *idx)
  {
    vectorGetElementIndices(idx);
    int index = getGlobalIndex();
    values[0] = gridpack::ComplexType(static_cast<double>(index),0.0);
    values[1] = gridpack::ComplexType(static_cast<double>(index+1),0.0);
  }

  void vectorSetElementValues(gridpack::ComplexType *values)
  {
    p_vec1 = values[0];
    p_vec2 = values[1];
  }

  void getValues(gridpack::ComplexType *values)
  {
    values[0] = p_vec1;
    values[1] = p_vec2;
  }

  void slabSize(int *isize, int *jsize)
  {
    *isize = 2;
    *jsize = NSLAB;
  }

  void slabSetRowIndex(int irow, int idx)
  {
    if (irow == 0) {
      p_slab_idx1 = idx;
    } else {
      p_slab_idx2 = idx;
    }
  }

  void slabGetRowIndices(int *idx)
  {
    idx[0] = p_slab_idx1;
    idx[1] = p_slab_idx2;
  }

  void slabGetValues(std::vector<gridpack::ComplexType*> &values, int *idx)
  {
    slabGetRowIndices(idx);
    int index = getGlobalIndex();
    int i;
    (values[0])[0] = gridpack::ComplexType(static_cast<double>(index),0.0);
    (values[1])[0] = gridpack::ComplexType(static_cast<double>(index+1),0.0);
    for (i=1; i<NSLAB; i++) {
      (values[0])[i] = (values[0])[i-1]+gridpack::ComplexType(1.0,0.0);
      (values[1])[i] = (values[1])[i-1]+gridpack::ComplexType(1.0,0.0);
    }
  }

  void slabSetValues(gridpack::ComplexType **values)
  {
    int k;
    for (k=0; k<NSLAB; k++) {
      (p_vals[0])[k] = real((values[0])[k]);
      (p_vals[1])[k] = real((values[1])[k]);
    }
  }

  double* getSlabArray(int idx)
  {
    return p_vals[idx];
  }

  gridpack::ComplexType p_val;
  gridpack::RealType p_rval;
  int p_row_idx;
  int p_col_idx;
  int p_vec_idx1;
  int p_vec_idx2;
  int p_slab_idx1;
  int p_slab_idx2;
  gridpack::ComplexType p_vec1, p_vec2;
  double **p_vals;
};

class TestBranch
  : public gridpack::component::BaseBranchComponent {
  public: 

  TestBranch(void) {
    p_vals = new double*[1];
    p_vals[0] = new double[NSLAB];
  }

  ~TestBranch(void) {
    delete [] p_vals[0];
    delete [] p_vals;
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

  int matrixNumRows() const {
    return 1;
  }

  int matrixNumCols() const {
    return 1;
  }

  void matrixSetRowIndex(int irow, int idx) {
    p_row_idx = idx;
  }

  void matrixSetColIndex(int icol, int idx) {
    p_col_idx = idx;
  }

  int matrixGetRowIndex(int idx) {
    return p_row_idx;
  }

  int matrixGetColIndex(int idx) {
    return p_col_idx;
  }

  int matrixNumValues() const {
    return 2;
  }

  void matrixGetValues(gridpack::ComplexType *values, int *rows, int *cols) {
    boost::shared_ptr<gridpack::component::BaseComponent> bus;
    int im, jm, idx, jdx1, jdx2;
    im = matrixGetRowIndex(0);
    bus = getBus1();
    jm = bus->matrixGetColIndex(0);
    values[0] = static_cast<gridpack::ComplexType>
      (dynamic_cast<gridpack::component::BaseBusComponent*>(bus.get()) ->getGlobalIndex());
    rows[0] = im;
    cols[0] = jm;
    bus = getBus2();
    jm = bus->matrixGetColIndex(0);
    values[1] = static_cast<gridpack::ComplexType>
      (dynamic_cast<gridpack::component::BaseBusComponent*>(bus.get())->getGlobalIndex());
    rows[1] = im;
    cols[1] = jm;
  }

  int vectorNumElements() const
  {
    return 1;
  }

  void vectorSetElementIndex(int ielem, int idx)
  {
    p_vec_idx = idx;
  }

  void vectorGetElementIndices(int *idx)
  {
    idx[0] = p_vec_idx;
  }

  void vectorGetElementValues(gridpack::ComplexType *values, int *idx)
  {
    vectorGetElementIndices(idx);
    int idx1 = getBus1GlobalIndex();
    int idx2 = getBus2GlobalIndex();
    values[0] = gridpack::ComplexType(static_cast<double>(idx1+idx2),0.0);
  }

  void vectorSetElementValues(gridpack::ComplexType *values)
  {
    p_vec_val = values[0];
  }

  void getValues(gridpack::ComplexType *values)
  {
    values[0] = p_vec_val;
  }

  void slabSize(int *isize, int *jsize)
  {
    *isize = 1;
    *jsize = NSLAB;
  }

  void slabSetRowIndex(int irow, int idx)
  {
    p_slab_idx = idx;
  }

  void slabGetRowIndices(int *idx)
  {
    idx[0] = p_slab_idx;
  }

  void slabGetValues(std::vector<gridpack::ComplexType*> &values, int *idx)
  {
    slabGetRowIndices(idx);
    int idx1 = getBus1GlobalIndex();
    int idx2 = getBus2GlobalIndex();
    int i;
    (values[0])[0] = gridpack::ComplexType(static_cast<double>(idx1+idx2),0.0);
    for (i=1; i<NSLAB; i++) {
      (values[0])[i] = (values[0])[i-1]+gridpack::ComplexType(1.0,0.0);
    }
  }

  void slabSetValues(gridpack::ComplexType **values)
  {
    int k;
    for (k=0; k<NSLAB; k++) {
      (p_vals[0])[k] = real((values[0])[k]);
    }
  }

  double* getSlabArray(int idx)
  {
    return p_vals[idx];
  }

  int p_row_idx;
  int p_col_idx;
  int p_vec_idx;
  int p_slab_idx;
  gridpack::ComplexType p_vec_val;
  double **p_vals;
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
  // Use bus_index array to keep track of local index of buses
  int *bus_index = new int[nx*ny];
  for (j=0; j<ny; j++) {
    iy = j + iaymin;
    for (i=0; i<nx; i++) {
      ix = i + iaxmin;
      n = iy*XDIM + ix;
      n = 2*n;  // Provide original index that is not equal to global index 
      // Add bus if new bus is attached to at least on local bus (bus is not
      // located at an interior corner of the grid)
      bool bus_added = false;
      if (ix == iaxmax && iy == iaymax) {
        if (iaxmax == ixmax || iaymax == iymax) {
          network->addBus(n);
          bus_added = true;
        }
      } else if (ix == iaxmin && iy == iaymax) {
        if (iaxmin == ixmin || iaymax == iymax) {
          network->addBus(n);
          bus_added = true;
        }
      } else if (ix == iaxmax && iy == iaymin) {
        if (iaxmax == ixmax || iaymin == iymin) {
          network->addBus(n);
          bus_added = true;
        }
      } else if (ix == iaxmin && iy == iaymin) {
        if (iaxmin == ixmin || iaymin == iymin) {
          network->addBus(n);
          bus_added = true;
        }
      } else {
        network->addBus(n);
        bus_added = true;
      }
      // Set active flag for network buses
      if (bus_added) {
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
        bus_index[j*nx+i] = ncnt;
        ncnt++;
      }
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
      n1 = bus_index[n1];
      n2 = ly*(iaxmax-iaxmin+1) + lx + 1;
      n2 = bus_index[n2];
      network->setLocalBusIndex1(ncnt,n1);
      network->setLocalBusIndex2(ncnt,n2);
      network->addBranchNeighbor(n1,ncnt);
      network->addBranchNeighbor(n2,ncnt);
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
      n1 = bus_index[n1];
      n2 = (ly+1)*(iaxmax-iaxmin+1) + lx;
      n2 = bus_index[n2];
      network->setLocalBusIndex1(ncnt,n1);
      network->setLocalBusIndex2(ncnt,n2);
      network->addBranchNeighbor(n1,ncnt);
      network->addBranchNeighbor(n2,ncnt);
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
  delete [] bus_index;

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
      if (network->getBus(i)->matrixDiagSize(&isize,&jsize)
          && isize > 0 && jsize > 0) {
        network->getBus(i)->getMatVecIndex(&idx);
        idx--;
        M->getElement(idx,idx,v);
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
    if (network->getBranch(i)->matrixForwardSize(&isize,&jsize)
        && isize > 0 && jsize > 0) {
      network->getBranch(i)->getMatVecIndices(&idx,&jdx);
      idx--;
      jdx--;
      if (idx >= rlo-1 && idx <= rhi-1) {
        M->getElement(idx,jdx,v);
        rv = real(v);
        if (rv != 1.0) {
          printf("p[%d] Forward matrix error i: %d j:%d v: %f\n",me,idx,jdx,rv);
          chk = 1;
        }
      }
      if (jdx >= rlo-1 && jdx <= rhi-1) {
        M->getElement(jdx,idx,v);
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
  boost::shared_ptr<gridpack::math::RealVector> rV = vMap.mapToRealVector();

  // Check to see if vector has correct values
  chk = 0;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      if (network->getBus(i)->vectorSize(&isize)) {
        network->getBus(i)->getMatVecIndex(&idx);
        idx--;
        V->getElement(idx,v);
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

  // Check to see if real vector has correct values
  chk = 0;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      if (network->getBus(i)->vectorSize(&isize)) {
        network->getBus(i)->getMatVecIndex(&idx);
        idx--;
        rV->getElement(idx,rv);
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
      printf("\nReal vector elements are ok\n");
    } else {
      printf("\nError found in real vector elements\n");
    }
  }

  if (me == 0) {
    printf("\nTesting mapToBus\n");
  }

  // Multiply values in vector by a factor of 2
  int lo, hi;
  V->localIndexRange(lo,hi);
  for (i=lo; i<hi; i++) {
    V->getElement(i,v);
    v *= 2.0;
    V->setElement(i,v);
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
          printf("p[%d] Bus error i: %d v: %f expected: %f\n",me,idx,rv,(double)(2*idx));
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
      printf("\nError found in bus values\n");
    }
  }

  if (me == 0) {
    printf("\nTesting general matrix interface\n");
  }
  chk = 0;
  gridpack::mapper::GenMatrixMap<TestNetwork> gMap(network); 
  boost::shared_ptr<gridpack::math::Matrix> G = gMap.mapToMatrix();
  int im, jm, nvals, icnt, jdx1, jdx2;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      im = network->getBus(i)->matrixGetRowIndex(0);
      nvals = network->getBus(i)->matrixNumValues();
      idx = network->getGlobalBusIndex(i);
      gridpack::ComplexType *values;
      values = new gridpack::ComplexType[nvals];
      icnt = 0;
      std::vector<int> nghbrs = network->getConnectedBranches(i);
      for (j=0; j<nghbrs.size(); j++) {
        if (network->getActiveBranch(nghbrs[j])) {
          network->getBranchEndpoints(nghbrs[j],&jdx1,&jdx2);
          if (i==jdx1) {
            jdx = network->getGlobalBusIndex(jdx2);
          } else {
            jdx = network->getGlobalBusIndex(jdx1);
          }
          jm = network->getBranch(nghbrs[j])->matrixGetColIndex(0);
          G->getElement(im,jm,v);
          if (v != static_cast<gridpack::ComplexType>(idx+jdx)) {
            chk = 1;
          }
        }
      }
      delete [] values;
      G->getElement(im,im,v);
      if (v != static_cast<gridpack::ComplexType>(idx)) {
         chk = 1;
      }
    }
  }

  GA_Igop(&chk,one,"+");
  if (me == 0) {
    if (chk == 0) {
      printf("\nGeneralized matrix values are ok\n");
    } else {
      printf("\nError found in generalized matrix value\n");
    }
  }

  if (me == 0) {
    printf("\nTesting generalized vector interface\n");
  }
  gridpack::mapper::GenVectorMap<TestNetwork> gvMap(network); 
  boost::shared_ptr<gridpack::math::Vector> GV = gvMap.mapToVector();

  // Check to see if vector has correct values
  chk = 0;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      idx = network->getBus(i)->getGlobalIndex();
      int indices[2];
      network->getBus(i)->vectorGetElementIndices(indices);
      GV->getElement(indices[0],v);
      if (idx != static_cast<int>(real(v))) {
        printf("Mismatch found (bus) expected: %d actual: %d\n",
            idx,static_cast<int>(real(v)));
        chk = 1;
      }
      GV->getElement(indices[1],v);
      if (idx+1 != static_cast<int>(real(v))) {
        printf("Mismatch found (bus) expected: %d actual: %d\n",
            idx+1,static_cast<int>(real(v)));
        chk = 1;
      }
    }
  }
  for (i=0; i<nbranch; i++) {
    if (network->getActiveBranch(i)) {
      int idx1 = network->getBranch(i)->getBus1GlobalIndex();
      int idx2 = network->getBranch(i)->getBus2GlobalIndex();
      int index;
      network->getBranch(i)->vectorGetElementIndices(&index);
      GV->getElement(index,v);
      if (idx1+idx2 != static_cast<int>(real(v))) {
        printf("Mismatch found (branch) expected: %d actual: %d\n",
            idx1+idx2,static_cast<int>(real(v)));
        chk = 1;
      }
    }
  }
  GA_Igop(&chk,one,"+");
  if (me == 0) {
    if (chk == 0) {
      printf("\nGeneralized vector elements are ok\n");
    } else {
      printf("\nError found in generalized vector elements\n");
    }
  }

  if (me == 0) {
    printf("\nTesting mapToNetwork for general vector\n");
  }

  // Push values back onto network
  gvMap.mapToNetwork(GV);

  // Check to see if network components have correct values
  chk = 0;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      idx = network->getBus(i)->getGlobalIndex();
      gridpack::ComplexType values[2];
      network->getBus(i)->getValues(values);
      if (idx != static_cast<int>(real(values[0]))) {
        printf("Mismatch found (bus) expected: %d actual: %d\n",
            idx,static_cast<int>(real(values[0])));
        chk = 1;
      }
      if (idx+1 != static_cast<int>(real(values[1]))) {
        printf("Mismatch found (bus) expected: %d actual: %d\n",
            idx+1,static_cast<int>(real(values[1])));
        chk = 1;
      }
    }
  }
  for (i=0; i<nbranch; i++) {
    if (network->getActiveBranch(i)) {
      int idx1 = network->getBranch(i)->getBus1GlobalIndex();
      int idx2 = network->getBranch(i)->getBus2GlobalIndex();
      gridpack::ComplexType value;
      int index;
      network->getBranch(i)->getValues(&value);
      if (idx1+idx2 != static_cast<int>(real(value))) {
        printf("Mismatch found (branch) expected: %d actual: %d\n",
            idx1+idx2,static_cast<int>(real(value)));
        chk = 1;
      }
    }
  }
  GA_Igop(&chk,one,"+");
  if (me == 0) {
    if (chk == 0) {
      printf("\nNetwork values are ok\n");
    } else {
      printf("\nError found in network value\n");
    }
  }

  if (me == 0) {
    printf("\nTesting generalized slab interface\n");
  }
  gridpack::mapper::GenSlabMap<TestNetwork> gsMap(network); 
  boost::shared_ptr<gridpack::math::Matrix> GS = gsMap.mapToMatrix();

  // Check to see if matrix has correct values
  chk = 0;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      idx = network->getBus(i)->getGlobalIndex();
      int indices[2];
      network->getBus(i)->slabGetRowIndices(indices);

      for (j=0; j<NSLAB; j++) {
        GS->getElement(indices[0],j,v);
        if (idx+j != static_cast<int>(real(v))) {
          printf("Mismatch found (bus) expected: %d actual: %d\n",
              idx+j,static_cast<int>(real(v)));
          chk = 1;
        }
        GS->getElement(indices[1],j,v);
        if (idx+j+1 != static_cast<int>(real(v))) {
          printf("Mismatch found (bus) expected: %d actual: %d\n",
              idx+j+1,static_cast<int>(real(v)));
          chk = 1;
        }
      }
    }
  }
  for (i=0; i<nbranch; i++) {
    if (network->getActiveBranch(i)) {
      int idx1 = network->getBranch(i)->getBus1GlobalIndex();
      int idx2 = network->getBranch(i)->getBus2GlobalIndex();
      int index;
      network->getBranch(i)->vectorGetElementIndices(&index);
      for(j=0; j<NSLAB; j++) {
        GS->getElement(index,j,v);
        if (idx1+idx2+j != static_cast<int>(real(v))) {
          printf("Mismatch found (branch) expected: %d actual: %d\n",
              idx1+idx2+j,static_cast<int>(real(v)));
          chk = 1;
        }
      }
    }
  }
  GA_Igop(&chk,one,"+");
  if (me == 0) {
    if (chk == 0) {
      printf("\nSlab matrix elements are ok\n");
    } else {
      printf("\nError found in slab matrix elements\n");
    }
  }

  // Test mapToNetwork function for a slab mapper
  if (me == 0) {
    printf("\nTesting mapToNetwork function for  slab interface\n");
  }
  gsMap.mapToNetwork(*GS);
  // Check to see if network has correct values
  chk = 0;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      idx = network->getBus(i)->getGlobalIndex();
      double *d0 = network->getBus(i)->getSlabArray(0);
      double *d1 = network->getBus(i)->getSlabArray(1);
      for (j=0; j<NSLAB; j++) {
        if (idx+j != static_cast<int>(d0[j])) {
          printf("Mismatch found (bus) for array 0 expected: %d actual: %d\n",
              idx+j,static_cast<int>(d0[j]));
          chk = 1;
        }
        if (idx+j+1 != static_cast<int>(d1[j])) {
          printf("Mismatch found (bus) for array 1 expected: %d actual: %d\n",
              idx+j+1,static_cast<int>(d1[j]));
          chk = 1;
        }
      }
    }
  }
  for (i=0; i<nbranch; i++) {
    if (network->getActiveBranch(i)) {
      int idx1 = network->getBranch(i)->getBus1GlobalIndex();
      int idx2 = network->getBranch(i)->getBus2GlobalIndex();
      double *d = network->getBranch(i)->getSlabArray(0);
      for(j=0; j<NSLAB; j++) {
        if (idx1+idx2+j != static_cast<int>(d[j])) {
          printf("Mismatch found (branch) for array expected: %d actual: %d\n",
              idx1+idx2+j,static_cast<int>(d[j]));
          chk = 1;
        }
      }
    }
  }
  GA_Igop(&chk,one,"+");
  if (me == 0) {
    if (chk == 0) {
      printf("\nSlab mapToNetwork ok\n");
    } else {
      printf("\nError found in array from matrixToNetwork\n");
    }
  }

}

int
main (int argc, char **argv) {

  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);

  // Initialize GA libraries
  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  // Initialize Math libraries
  gridpack::math::Initialize(&argc,&argv);
  gridpack::parallel::Communicator comm;
  MPI_Comm world = static_cast<MPI_Comm>(comm);

  int me;
  ierr = MPI_Comm_rank(world, &me);
  int nprocs;
  ierr = MPI_Comm_size(world, &nprocs);
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
  return ierr;
}
