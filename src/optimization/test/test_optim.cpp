/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#include <vector>
#include <iostream>
#include <cstring>
#include <stdio.h>


#include "mpi.h"
#include <macdecls.h>
#include "gridpack/utilities/complex.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/optimization/optimization.hpp"
#include <ilcplex/ilocplex.h>
#include <stdlib.h>

#define XDIM 1
#define YDIM 1
#define NSLAB 20

class TestBus
  : public gridpack::component::BaseBusComponent {
  public: 
    /**
     *  Simple constructor
     */
    TestBus(void) {
//      p_numUnits = 0;
//      p_numHorizons = 0;
    }
    /**
     *  Simple destructor
     */
    ~TestBus(void) {
    }
    /**
     * Load values stored in DataCollection object into HWBus object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       bus that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data)
    {
//       data->getValue(BUS_,&p_original_idx);

    }

    double objectiveFunction(void) {
      double obj;
      obj = 2.0;
      return obj;
    }

    /**
     * Return the size of the generator on the network.
     */
//    int netGenSize(int *isize) {
//      *isize = p_numUnit;
//      GA_Igop(isize,one,"+");
//    }
    /**
     * Return the size of the buffer used in data exchanges on the network.
     * For this problem, the number of plant units need to be exchanged
     * @return size of buffer
     */
//  void getXCBufSize(void);
    /**
     * Assign pointers for plant units
     */
//  void setXCBuf(void *buf);

    /**
     * Set values of constant cost parameters. These can then be used in subsequent
     * calculations of cost function
     */
    void setCostConst(void);
    /**
     * Get values of constant cost parameters. These can then be used in subsequent
     * calculations
     */
    void getCostConst(void);

    /**
     * Set values of linear cost parameters. These can then be used in subsequent
     * calculations of cost function
     */
    void setLinearConst(void);
    /**
     * Get values of linear cost parameters. These can then be used in subsequent
     * calculations
     */
    void getLinearConst(void);
    /**
     * Set values of quadratic cost parameters. These can then be used in subsequent
     * calculations of cost function
     */
    void setQuadConst(void);
    /**
     * Get values of quadratic cost parameters. These can then be used in subsequent
     * calculations
     */
    void getQuadConst(void);

//    gridpack::ComplexType p_val;
//    int p_row_idx;
//    int p_col_idx;
//    int p_vec_idx1;
//    int p_vec_idx2;
//    int p_slab_idx1;
//    int p_slab_idx2;
//    gridpack::ComplexType p_vec1, p_vec2;
  private:
//    int p_numUnits;

//  friend class boost::serialization::access;

//  template<class Archive>
//  void serialize(Archive & ar, const unsigned int version)
//  {
//    ar & boost::serialization::base_object<gridpack::component::BaseBusComponent>(*this)
//      & p_original_idx;
//  }

};

class TestBranch
  : public gridpack::component::BaseBranchComponent {
  public: 

    TestBranch(void) {
    }

    ~TestBranch(void) {
    }

    int p_row_idx;
    int p_col_idx;
    int p_vec_idx;
    int p_slab_idx;
    gridpack::ComplexType p_vec_val;
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

int p_numUnits;
std::vector<double> p_minPower;
std::vector<double> p_maxPower;
std::vector<int> p_minUpTime;
std::vector<int> p_minDownTime;
std::vector<double> p_costConst;
std::vector<double> p_costLinear;
std::vector<double> p_costQuad;

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
// set up generator info on each bus
    int t_numUnits;
    int t_numHorizons;
    std::vector<double> t_demand;
    std::vector<double> t_minPower;
    std::vector<double> t_maxPower;
    std::vector<int> t_minUptime;
    std::vector<int> t_minDowntime;
    std::vector<double> t_costConst;
    std::vector<double> t_costLinear;
    std::vector<double> t_costQuad;
  gridpack::component::DataCollection *data;
// 1 bus
//  for (i=0; i<nbus; i++) {
        int l_idx = 0;
        data = dynamic_cast<gridpack::component::DataCollection*>
          (network->getBusData(l_idx).get());
        t_numUnits = 3;
    if (!data->setValue("GENERATOR_NUMBERS",t_numUnits)) {
        data->addValue("GENERATOR_NUMBERS",t_numUnits);
    }
        t_numHorizons = 3;
        data->setValue("GENERATOR_HORIZONS",t_numHorizons);
        t_demand.push_back(150);
        t_demand.push_back(300);
        t_demand.push_back(200);
        t_minUptime.push_back(2);
        t_minUptime.push_back(1);
        t_minUptime.push_back(1);
        t_minDowntime.push_back(3);
        t_minDowntime.push_back(1);
        t_minDowntime.push_back(2);
        t_costConst.push_back(510.00);
        t_costConst.push_back(310.0);
        t_costConst.push_back ( 78.0);
        t_costLinear.push_back(7.9);
        t_costLinear.push_back(7.85);
        t_costLinear.push_back(9.56);
        t_costQuad.push_back(0.00172);
        t_costQuad.push_back(0.00194);
        t_costQuad.push_back(0.00694);
        t_minPower.push_back(150.0);
        t_minPower.push_back(50);
        t_minPower.push_back(10.0);
        t_maxPower.push_back(250.0);
        t_maxPower.push_back(100.0);
        t_maxPower.push_back( 50.0);
  double rval;
  int ival;
// The following should be read from a network
  for (i=0; i<t_numUnits; i++) {
    rval = t_demand[i];
    if (!data->setValue("GENERATOR_DEMAND",rval,i)) {
      data->addValue("GENERATOR_DEMAND",rval,i);
    }
    ival = t_minUptime[i];
    if (!data->setValue("GENERATOR_MIN_UPTIME",ival,i)) {
      data->addValue("GENERATOR_MIN_UPTIME",ival,i);
    }
    ival = t_minDowntime[i];
    if (!data->setValue("GENERATOR_MIN_DNTIME",ival,i)) {
      data->addValue("GENERATOR_MIN_DNTIME",ival,i);
    }
    rval = t_minPower[i];
    if (!data->setValue("GENERATOR_MIN_POWER",rval,i)) {
      data->addValue("GENERATOR_MIN_POWER",rval,i);
    }
    rval = t_maxPower[i];
    if (!data->setValue("GENERATOR_MAX_POWER",rval,i)) {
      data->addValue("GENERATOR_MAX_POWER",rval,i);
    }
    rval = t_costConst[i];
    if (!data->setValue("GENERATOR_CONST_COST",rval,i)) {
      data->addValue("GENERATOR_CONST_COST",rval,i);
    }
    rval = t_costLinear[i];
    if (!data->setValue("GENERATOR_LINEAR_COST",rval,i)) {
      data->addValue("GENERATOR_LINEAR_COST",rval,i);
    }
    rval = t_costQuad[i];
    if (!data->setValue("GENERATOR_QUAD_COST",rval,i)) {
      data->addValue("GENERATOR_QUAD_COST",rval,i);
    }
  }


  // Set up some remaining properties of network and network components so that
  // matrix-vector interface is ready to go
  gridpack::factory::BaseFactory<TestNetwork> factory(network);
  factory.setComponents();
  factory.load();

  if (me == 0) {
    printf("\nTesting Optimizer\n");
  }

  
  // Check to see if objective function has correct values
  int one = 1;
  int chk = 0;
  int nbus = network->numBuses();
  int nbranch = network->numBranches();
  int idx,jdx,isize,jsize;
  double obj,objsum;
  double rv;
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
        obj = network->getBus(i)->objectiveFunction();
        if (obj != 2.0) {
          printf("p[%d] Objection function error i: %d j:%d v: %f\n",me,i,i,obj);
          chk = 1;
        }
    }
  }
  GA_Igop(&chk,one,"+");
  if (me == 0) {
    if (chk == 0) {
      printf("\nOjectvie functions are ok\n");
    } else {
      printf("\nError found in objective functions\n");
    }
  }
  gridpack::optimization::NetworkOptimizer<TestNetwork> optim(network);
// Get data from the network
  optim.loadBusData(); 
  p_numUnits = optim.numUnits;
  p_minPower = optim.minPower;
  p_maxPower = optim.maxPower;
  p_minUpTime = optim.minUpTime;
  p_minDownTime = optim.minDownTime;
  p_costConst = optim.costConst;
  p_costLinear = optim.costLinear;
  p_costQuad = optim.costQuad;
//  printf ("num of unit --, %f \n",optim.minPower[0]);

//Get the objective function on the network
  objsum = optim.objectiveFunction();  
  if (me == 0) {
      printf("\nNetwork ojectvie function is: %f\n\n",objsum);
  }
 
}


//A test unit commitment problem
typedef IloArray<IloIntVarArray> IntArray2;
typedef IloArray<IloNumVarArray> NumArray2;

void run_unit_commitment() {
    
//  gridpack::optimization::Optimizer<TestNetwork> optim(network);
// Get data from the network
//  optim.loadBusData(); 
    double rval;
    int ival;
    
    IloEnv env;
    //
    // Data
    //
    const IloInt numUnits = p_numUnits;
    const IloInt numHorizons = 3;
    // create minimum power array and append values to the array
    //IloNumArray minPower(env,numUnits,0);
    //minPower = p_minPower;
    IloNumArray minPower(env);
    IloNumArray maxPower(env);
    IloIntArray minUpTime(env);
    IloIntArray minDownTime(env);
    IloNumArray costConst(env);
    IloNumArray costLinear(env);
    IloNumArray costQuad(env);
    for (int i = 0; i < numUnits; i++) {
      rval = p_minPower[i];
      minPower.add(rval);
    }
    for (int i = 0; i < numUnits; i++) {
      rval = p_maxPower[i];
      maxPower.add(rval);
    }
    for (int i = 0; i < numUnits; i++) {
      ival = p_minUpTime[i];
      minUpTime.add(ival);
    }
    for (int i = 0; i < numUnits; i++) {
      ival = p_minDownTime[i];
      minDownTime.add(ival);
    }
    for (int i = 0; i < numUnits; i++) {
      rval = p_costConst[i];
      costConst.add(rval);
    }
    for (int i = 0; i < numUnits; i++) {
      rval = p_costLinear[i];
      costLinear.add(rval);
    }
    for (int i = 0; i < numUnits; i++) {
      rval = p_costQuad[i];
      costQuad.add(rval);
    }
    // marginal cost and startup cost
    IloNumArray margCost(env, numUnits, 10.0, 12.0, 20.0);
    IloNumArray startupCost(env, numUnits, 1000.0, 600.0, 100.0);
    // Demand at each horizon
    IloIntArray demand(env, numHorizons, 150, 300, 200);
    // Parameters for cost function
    // Variable arrays
    IntArray2 onOff;
    NumArray2 powerProduced;
    onOff = IntArray2(env,numHorizons);
    for (IloInt p = 0; p < numHorizons; p++) {
      onOff[p] = IloIntVarArray(env, numUnits,0,1);
    }
    powerProduced = NumArray2(env,numHorizons);
    for (IloInt p = 0; p < numHorizons; p++) {
      powerProduced[p] = IloNumVarArray(env, 0,maxPower,ILOFLOAT);
    }
//  Objective
    IloModel ucmdl(env);
    IloExpr obj(env);
    for (IloInt p = 0; p < numHorizons; p++) {
      for (IloInt i = 0; i < numUnits; i++) {
         obj += costConst[i]*onOff[p][i] +
               costLinear[i]*powerProduced[p][i] 
              + costQuad[i]*powerProduced[p][i]*powerProduced[p][i];
      }
    }
    
    ucmdl.add(IloMinimize(env,obj));
    obj.end();


//
//  Constraints
//  
//  Initial state, treat as constraint
    ucmdl.add(onOff[0][0] == 1);
    ucmdl.add(powerProduced[0][0] == 150);
    IloInt upDnPeriod;
    for (IloInt p = 1; p < numHorizons; p++) {
      IloExpr expr3(env);
      for (IloInt i = 0; i < numUnits; i++) {
         IloExpr expr1(env);
         IloExpr expr2(env);
         expr1 = powerProduced[p][i] - 10000*onOff[p][i];
         ucmdl.add( expr1 <= 0);
         expr2 = powerProduced[p][i] - minPower[i]*onOff[p][i];
         ucmdl.add( expr2 >= 0);
         expr1.end();
         expr2.end();
// minium up and down time
// on at horizon p
         IloExpr upDnIndicator(env);
         upDnIndicator = onOff[p][i] - onOff[p-1][i]; 
         upDnPeriod = std::min(numHorizons, (p+minUpTime[i]));
         for (IloInt j = p; j < upDnPeriod; j++) {
           ucmdl.add( upDnIndicator - 10000*onOff[j][i] <= 0);
         }
         upDnIndicator.end();
// off at horizon p
         upDnIndicator = 1 - (onOff[p-1][i] - onOff[p][i]); 
         upDnPeriod = std::min(numHorizons, (p+minDownTime[i]));
         for (IloInt j = p; j < upDnPeriod; j++) {
           ucmdl.add( upDnIndicator - 10000*onOff[j][i] <= 0);
         }
         upDnIndicator.end();
      }
      expr3 = IloSum(powerProduced[p]);
      ucmdl.add( expr3 == demand[p]);
      expr3.end();
    }
// Solve model
    IloCplex cplex(env);
// setup parallel mode and number of threads
    cplex.setParam(IloCplex::ParallelMode, 1);
    cplex.setParam(IloCplex::Threads, 2);
    cplex.extract(ucmdl);
    cplex.exportModel("test_uc.lp");
    std::ofstream logfile("cplex.log");
    cplex.setOut(logfile);
    cplex.setWarning(logfile);
    cplex.solve();
    cplex.out() << "solution status = " << cplex.getStatus() << std::endl;
    cplex.out() << "cost   = " << cplex.getObjValue() << std::endl;
    for (IloInt p = 0; p < numHorizons; p++) {
      for (IloInt i = 0; i < numUnits; i++) {
          env.out() << "At time " << p << " Power produced by unit " << i << " " <<
          cplex.getValue(onOff[p][i]) << "  " <<
          cplex.getValue(powerProduced[p][i]) << std::endl;
      }
    }

  
  env.end();
}

int
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
    printf("Testing Optimizer Module\n");
    printf("\nTest Network is %d X %d\n",XDIM,YDIM);
  }

  run(me, nprocs);
  run_unit_commitment();

  GA_Terminate();

  // Clean up MPI libraries
  ierr = MPI_Finalize();
}
