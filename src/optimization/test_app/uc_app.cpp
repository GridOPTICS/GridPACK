/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   uc_app.cpp
 * @author 
 * @date   2015-11-03 13:40:11 d3g096
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

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp> // needed of is_any_of()
#include "uc_optimizer.hpp"
#include "uc_factory.hpp"
#include "uc_app.hpp"
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/parser/hash_distr.hpp"
#include "mpi.h"
#include <stdlib.h>
#include <ga.h>
#include "gridpack/network/base_network.hpp"
#include "gridpack/optimization/optimizer.hpp"
#include <boost/format.hpp>
typedef gridpack::unit_commitment::UCBus::uc_ts_data uc_ts_data;
int Horizons;

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
 * Get time series data for loads and reserves and move them to individual
 * buses
 * @param filename name of file containing load and reserves time series
 * data
 */
void gridpack::unit_commitment::UCApp::getLoadsAndReserves(const char* filename)
{
  int me = p_network->communicator().rank();
  std::ifstream input;
  int nvals = 0;
  int i;
  gridpack::utility::StringUtils util;
  std::string line;
  std::vector<std::string>  split_line;
  if (me == 0) {
    input.open(filename);
    if (!input.is_open()) {
      char buf[512];
      sprintf(buf,"Failed to open file with load and reserve time series data: %s\n\n",
          filename);
      throw gridpack::Exception(buf);
    }
    std::getline(input,line);
    boost::split(split_line, line, boost::is_any_of("\t "),boost::token_compress_on );
    nvals = split_line.size();
    nvals = (nvals-1)/2;
  }
  p_network->communicator().sum(&nvals,1);
  std::vector<uc_ts_data*> series; 
  std::vector<int> bus_ids;
  // read in times series values on each of the buses
  if (me == 0) {
    while (std::getline(input, line)) {
      uc_ts_data *data = new uc_ts_data[nvals];
      boost::split(split_line, line, boost::is_any_of("\t "),boost::token_compress_on);
      bus_ids.push_back(atoi(split_line[0].c_str()));
      for (i=0; i<nvals; i++) {
        data[i].load = atof(split_line[2*i+1].c_str());
        data[i].reserve = atof(split_line[2*i+2].c_str());
      }
      series.push_back(data);
    }
    input.close();
  }
  // distribute data from process 0 to the processes that have the corresponding
  // buses
  gridpack::hash_distr::HashDistribution<UCNetwork,uc_ts_data,uc_ts_data>
    hash(p_network);
  int ndata;
    ndata = nvals;
  Horizons = ndata;
  hash.distributeBusValues(bus_ids, series, nvals);
  nvals = bus_ids.size();
// number of bus with loads input
  std::vector<double> loads;
  std::vector<double> reserve;
  for (i=0; i<nvals; i++) {
    int l_idx = bus_ids[i];
    gridpack::unit_commitment::UCBus *bus;
    bus = dynamic_cast<gridpack::unit_commitment::UCBus*>(
        p_network->getBus(l_idx).get());
    bus->setTimeSeries(series[i], ndata);
  }
}

/**
 * Execute application
 * @param argc number of arguments
 * @param argv list of character strings
 */

void gridpack::unit_commitment::UCApp::execute(int argc, char** argv)
{
  // load input file
  gridpack::parallel::Communicator world;
  p_network.reset(new UCNetwork(world));

  // read configuration file
  std::string filename = "uc_test.raw";

  // Read in external PTI file with network configuration
  gridpack::parser::PTI23_parser<UCNetwork> parser(p_network);
  parser.parse(filename.c_str());

  // partition network
  p_network->partition();

  gridpack::unit_commitment::UCApp::getLoadsAndReserves("loads.csv");

  int me;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

    int p_nBuses = p_network->numBuses();
    int ngen;

    gridpack::component::DataCollection *data;
    gridpack::unit_commitment::UCBus *bus;
    
    double *uc_load;
    double *uc_reserve;
    uc_load = new double[Horizons] ();
    uc_reserve = new double[Horizons] ();

    for (int i_bus = 0; i_bus < p_nBuses; i_bus++) {
      if(p_network->getActiveBus(i_bus)) {
         data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(i_bus).get());
           bus = dynamic_cast<gridpack::unit_commitment::UCBus*>(
            p_network->getBus(i_bus).get());
// calcuate total load and reserve at each period
           if( (bus->p_load).size() > 0) {
             for (int i = 0; i < Horizons; i++) {
               uc_load[i] += bus->p_load[i];
               uc_reserve[i] += bus->p_reserve[i];
             }    
           }
//        }
      }
    }

    p_network->communicator().sum(uc_load,Horizons);
    p_network->communicator().sum(uc_reserve,Horizons);
  // load uc parameters
  filename = "gen.uc";
  parser.externalParse(filename.c_str());


  // create factory
  gridpack::unit_commitment::UCFactory factory(p_network);
  factory.load();

  // create optimization object
  gridpack::unit_commitment::UCoptimizer optim(p_network);
// get uc parameters
  optim.getUCparam();
// get demands and reserve for optimizer
  optim.getLoadsInfo(Horizons,uc_load,uc_reserve);
// get a vector of optimization variables
  typedef boost::shared_ptr<gridpack::optimization::Variable> VarPtr;
  typedef boost::shared_ptr<gridpack::optimization::Expression> ExpPtr;
  typedef boost::shared_ptr<gridpack::optimization::Constraint> ConstPtr;

  gridpack::optimization::Optimizer opt(world);

//return list of variables 
  std::vector<VarPtr> p_vlist;
  p_vlist.clear();
  p_vlist = optim.getVariables();
  for (std::vector<VarPtr>::iterator i = p_vlist.begin();
       i != p_vlist.end(); ++i) {
    opt.addVariable(*i);
  }

//return expression representing contribution to objective function.

//return list of constraints
  std::vector<ConstPtr> locConstraint;
  locConstraint = optim.getLocalConstraints();
//
  for (std::vector<ConstPtr>::iterator ic = locConstraint.begin();
       ic != locConstraint.end(); ++ic) {
    opt.addConstraint(*ic);
  }
//get the global constraints
  std::vector<ExpPtr> globalConstraint;
  globalConstraint = optim.getGlobalConstraint("Energy Balance");
  ExpPtr pempty;
  boost::format fmt("PBalance%i");
  ConstPtr pc;
  for (int i=1; i<Horizons; i++) {
    pc = pempty == uc_load[i] ;
    std::string cname = boost::str(fmt % i);
    pc->name(cname);
    opt.createGlobalConstraint(pc->name(), pc);
   }
  ExpPtr rempty;
  boost::format rfmt("RBalance%i");
  for (int i=1; i<Horizons; i++) {
    pc = pempty >= uc_reserve[i] ;
    std::string cname = boost::str(rfmt % i);
    pc->name(cname);
    opt.createGlobalConstraint(pc->name(), pc);
   }

//    opt.addToGlobalConstraint("PBalance",globalConstraint);
  int pcnt = 1;
  int rcnt = 1;
  int icnt = 1;
  for (std::vector<ExpPtr>::iterator ic = globalConstraint.begin();
       ic != globalConstraint.end(); ++ic) {
     if(icnt % 2 == 1) {
       std::string cname = boost::str(fmt % pcnt);
       opt.addToGlobalConstraint(cname, *ic);
       pcnt++;
     }
     else{
       std::string cname = boost::str(rfmt % rcnt);
       opt.addToGlobalConstraint(cname, *ic);
       rcnt++;
     }
     icnt++;
  }

  ExpPtr objFunc;
  objFunc = optim.getObjectiveFunction();
  opt.addToObjective(objFunc);
  opt.minimize();
} 
