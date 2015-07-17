/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   uc_app.cpp
 * @author 
 * @date   
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <vector>
#include <cstring>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>

#include "uc_optimizer.hpp"
#include "uc_factory.hpp"
#include "uc_app.hpp"
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "mpi.h"
#include <gurobi_c++.h>
#include <stdlib.h>
#include <deque>
using namespace std;



// Calling program for unit commitment optimization application

/**
 * Basic constructor
 */
gridpack::unit_commitment::UCApp::UCApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::unit_commitment::UCApp::~UCApp(void)
{
}

/**
 * Execute application
 * @param argc number of arguments
 * @param argv list of character strings
 */
//typedef IloArray<IloIntVarArray> IntArray2;
//typedef IloArray<IloNumVarArray> NumArray2;

void gridpack::unit_commitment::UCApp::execute(int argc, char** argv)
{
  // load input file
  gridpack::parallel::Communicator world;
  boost::shared_ptr<UCNetwork> network(new UCNetwork(world));

  // read configuration file
  std::string filename = "uc_test.raw";

  // Read in external PTI file with network configuration
  gridpack::parser::PTI23_parser<UCNetwork> parser(network);
  parser.parse(filename.c_str());

  // partition network
  network->partition();

  // load uc parameters
  filename = "gen.uc";
//  if (filename.size() > 0) parser.parse(filename.c_str());
  parser.parse(filename.c_str());
#if 0
  parser.externalParse(filename.c_str());
#endif


  // create factory
  gridpack::unit_commitment::UCFactory factory(network);
  factory.load();

  // create optimization object
  gridpack::unit_commitment::UCoptimizer optim(network);
  // get uc parameters
  optim.getUCparam();

  // prepare for optimization
  double rval;
  int ival;
  int me = MPI::COMM_WORLD.Get_rank();

  // create lp env
  GRBEnv *env = 0;
  env = new GRBEnv();
  if(me == 0) {
  //
  // read demands and reserve from an input file
  //
    std::vector<std::vector<double> > loads;
    std::ifstream fin("loads.txt");
    std::string line;
    while (std::getline(fin, line)) {
      std::vector<double> lineData;           // create a new row
      double val;
      std::istringstream lineStream(line); 
      while (lineStream >> val) {          // for each value in line
        lineData.push_back(val);           // add to the current row
      }
      loads.push_back(lineData);         // add row to loads
    }
    int size = loads.size();
    const int numHorizons = size;
//    const IloInt numHorizons = 2;
    const int numUnits = optim.totalGen;
//
// create array on each process to store results from optimization
    int arr_size = numUnits*numHorizons;
//   
    double *minPower;
    double *maxPower;
    double *minUpTime;
    double *minDownTime;
    double *costConst;
    double *costLinear;
    double *costQuad;
    double *iniLevel;
    double *shutCap;
    double *startCap;
    double *startUp;
    double *rampUp;
    double *rampDown;
    double *initPeriod;
    double *demand;
    double *reserve;
    double val;
    minPower = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      minPower[i] = optim.uc_minPower[i];
    }
    maxPower = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      maxPower[i] = optim.uc_maxPower[i];
    }
    minUpTime = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      minUpTime[i] = optim.uc_minUpTime[i];
    }
    minDownTime = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      minDownTime[i] = optim.uc_minDownTime[i];
    }
    costConst = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      costConst[i] = optim.uc_costConst[i];
    }
    costLinear = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      costLinear[i] = optim.uc_costLinear[i];
    }
    costQuad = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      costQuad[i] = optim.uc_costQuad[i];
    }
    iniLevel = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      iniLevel[i] = optim.uc_iniLevel[i];
    }
    startUp = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      startUp[i] = optim.uc_startUp[i];
    }
    startCap = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      startCap[i] = optim.uc_startCap[i];
    }
    shutCap = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      shutCap[i] = optim.uc_shutCap[i];
    }
    rampUp = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      rampUp[i] = optim.uc_rampUp[i];
    }
    rampDown = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      rampDown[i] = optim.uc_rampDown[i];
    }
    initPeriod = new double[numUnits] ();
    for (int i = 0; i < numUnits; i++) {
      initPeriod[i] = optim.uc_initPeriod[i];
    }
    demand = new double[numHorizons] ();
    for (int i = 0; i < numHorizons; i++) {
      demand[i] = loads[i][0];
    }
    reserve = new double[numHorizons] ();
    for (int i = 0; i < numHorizons; i++) {
      reserve[i] = loads[i][1];
    }
//
    GRBModel ucmdl = GRBModel(*env);
    GRBVar** onOff = 0;
    onOff = new GRBVar*[numHorizons];
    for (int p = 0; p < numHorizons; p++) {
      onOff[p] = new GRBVar[numUnits];
      for (int i = 0; i < numUnits; i++) {
        onOff[p][i] = ucmdl.addVar(0.0, 1.0, 0.0, GRB_BINARY);
      }
    }
// Update model to integrate new variables
    ucmdl.update();
//
    GRBVar** start_Up = 0;
    start_Up = new GRBVar*[numHorizons];
    for (int p = 0; p < numHorizons; p++) {
      start_Up[p] = new GRBVar[numUnits];
      for (int i = 0; i < numUnits; i++) {
        start_Up[p][i] = ucmdl.addVar(0.0, 1.0, 0.0, GRB_BINARY);
      }
    }
    ucmdl.update();
//
    GRBVar** shutDown = 0;
    shutDown = new GRBVar*[numHorizons];
    for (int p = 0; p < numHorizons; p++) {
      shutDown[p] = new GRBVar[numUnits];
      for (int i = 0; i < numUnits; i++) {
        shutDown[p][i] = ucmdl.addVar(0.0, 1.0, 0.0, GRB_BINARY);
      }
    }
    ucmdl.update();
//
    GRBVar** powerProduced = 0;
    powerProduced = new GRBVar*[numHorizons];
    for (int p = 0; p < numHorizons; p++) {
      powerProduced[p] = new GRBVar[numUnits];
      for (int i = 0; i < numUnits; i++) {
        powerProduced[p][i] = ucmdl.addVar(0.0, maxPower[i], 0.0, GRB_CONTINUOUS);
      }
    }
    ucmdl.update();
//
    GRBVar** powerReserved = 0;
    powerReserved = new GRBVar*[numHorizons];
    for (int p = 0; p < numHorizons; p++) {
      powerReserved[p] = ucmdl.addVars(numUnits);
      for (int i = 0; i < numUnits; i++) {
        powerReserved[p][i] = ucmdl.addVar(0.0, maxPower[i], 0.0, GRB_CONTINUOUS);
      }
    }
    ucmdl.update();
//  Objective
    GRBQuadExpr obj = 0;
    for (int p = 1; p < numHorizons; p++) {
      for (int i = 0; i < numUnits; i++) {
         obj += costConst[i]*onOff[p][i]
              + startUp[i]*start_Up[p][i]
              + costLinear[i]*powerProduced[p][i]
              + costQuad[i]*powerProduced[p][i]*powerProduced[p][i];
      }
    }
    ucmdl.setObjective(obj);
// The objective is to minimize total costs
//    ucmdl.set(GRB_IntAttr_ModelSense, 1);
    
//
//  Constraints
//  
//  Initial state, treat as constraint
    for(int i=0; i< numUnits; i++) {
      ucmdl.addConstr(onOff[0][i] == 1);
      ucmdl.addConstr(start_Up[0][i] == 0);
      ucmdl.addConstr(shutDown[0][i] == 0);
      ucmdl.addConstr(powerProduced[0][i] == iniLevel[i]);
    }
    int upDnPeriod;
    for (int p = 1; p < numHorizons; p++) {
      GRBLinExpr expr3;
      for (int i = 0; i < numUnits; i++) {
         GRBLinExpr expr1;
         GRBLinExpr expr2;
         expr1 = powerProduced[p][i] - 10000*onOff[p][i];
         ucmdl.addConstr( expr1 <= 0);
         expr2 = powerProduced[p][i] - minPower[i]*onOff[p][i];
         ucmdl.addConstr( expr2 >= 0);
// ramp up constraint
         expr1 = powerProduced[p][i]+powerReserved[p][i]-powerProduced[p-1][i];
         ucmdl.addConstr( expr1 <= rampUp[i]);
// ramp down constraint
         expr2 = powerProduced[p-1][i]-powerProduced[p][i];
         ucmdl.addConstr( expr2 <= rampDown[i]);
// minium up and down time
// on at horizon p
         GRBLinExpr upDnIndicator;
         upDnIndicator = onOff[p][i] - onOff[p-1][i];
         if(p == 1) {
          int initP = (int)(initPeriod[i]);
          GRBLinExpr upDnIndicator0;
//
          upDnIndicator0 = onOff[p-1][i]-onOff[p][i];
          upDnPeriod = std::min(numHorizons, (p+(int)(minUpTime[i]+0.5)-initP));
          for (int j = p; j < upDnPeriod; j++) {
             ucmdl.addConstr( upDnIndicator0 - 10000*onOff[j][i] <= 0);
          }
         } else{
           upDnPeriod = std::min(numHorizons, (p+(int)(minUpTime[i]+0.5)));
          for (int j = p; j < upDnPeriod; j++) {
           ucmdl.addConstr( upDnIndicator - 10000*onOff[j][i] <= 0);
          }
         }
// start up, previous off
         ucmdl.addConstr( upDnIndicator - 10000*start_Up[p][i] <= 0);
// off at horizon p
         upDnIndicator = 1 - (onOff[p-1][i] - onOff[p][i]);
         upDnPeriod = std::min(numHorizons, (p+(int)(minDownTime[i]+0.5)));
         for (int j = p; j < upDnPeriod; j++) {
           ucmdl.addConstr( upDnIndicator - 10000*onOff[j][i] <= 0);
         }
// shut down, previous on
         upDnIndicator = onOff[p-1][i] - onOff[p][i];
         ucmdl.addConstr( upDnIndicator - 10000*shutDown[p][i] <= 0);
// generation limits
// startup at horizon p
         expr1 = powerProduced[p][i]-minPower[i]+powerReserved[p][i];
         expr2 = (maxPower[i]-minPower[i])*onOff[p][i]
               -(maxPower[i]-startCap[i])*start_Up[p][i];
         ucmdl.addConstr( expr1 <= expr2);
// shutdown at horizon p 
         expr1 = powerProduced[p-1][i]-minPower[i]+powerReserved[p-1][i];
         expr2 = (maxPower[i]-minPower[i])*onOff[p-1][i]
               -(maxPower[i]-shutCap[i])*shutDown[p][i];
         ucmdl.addConstr( expr1 <= expr2);
      }
      expr3 = 0.0;
      for (int i = 0; i < numUnits; i++) {
        expr3 += powerProduced[p][i];
      }
      ucmdl.addConstr( expr3 == demand[p]);

      expr3 = 0.0;
      for (int i = 0; i < numUnits; i++) {
        expr3 += powerReserved[p][i];
      }
      ucmdl.addConstr( expr3 >= reserve[p]);
    }
    ucmdl.update();
//Outputmodel to a file
    ucmdl.write("test_gurobi.lp");

// Solve model
    ucmdl.optimize();

// Extract solution
    int status = ucmdl.get(GRB_IntAttr_Status);
    if ((status == GRB_INF_OR_UNBD) || (status == GRB_INFEASIBLE) ||
      (status == GRB_UNBOUNDED))
    {
      cout << "The model cannot be solved " <<
      "because it is infeasible or unbounded" << endl;
      return;
    }

    if (status != GRB_OPTIMAL)
    {
      cout << "Optimization was stopped with status " << status << endl;
      return;
    }

// write solution    
    int busid;
    for (int i = 0; i < numUnits; i++) {
      busid = optim.busID[i];
      for (int p = 0; p < numHorizons; p++) {
          cout << "At time " << p << " Power produced by unit " << i << " " << "on bus  " << busid << " " <<
          onOff[p][i].get(GRB_DoubleAttr_X) << "  " <<
          start_Up[p][i].get(GRB_DoubleAttr_X) << "  " <<
          shutDown[p][i].get(GRB_DoubleAttr_X) << "  " <<
          powerReserved[p][i].get(GRB_DoubleAttr_X) << "  " <<
          powerProduced[p][i].get(GRB_DoubleAttr_X) << std::endl;
      }
    }

  for (int p=0; p < numHorizons; p++) {
    delete[] onOff[p];
    delete[] start_Up[p];
    delete[] shutDown[p];
    delete[] powerReserved[p];
    delete[] powerProduced[p];
  }
  delete[] onOff;
  delete[] start_Up;
  delete[] shutDown;
  delete[] powerReserved;
  delete[] powerProduced;
  }
  delete env;
} 
