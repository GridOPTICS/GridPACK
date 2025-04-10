/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_app_module.cpp
 * @author Yousu Chen, Bruce Palmer
 * @date   2014-09-18 12:27:18 d3g096
 * Last updated: 8/5/2014 
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#include "gridpack/configuration/configuration.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/parser/PTI33_parser.hpp"
#include "gridpack/parser/PTI34_parser.hpp"
#include "gridpack/parser/PTI35_parser.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/gen_matrix_map.hpp"
#include "gridpack/mapper/gen_vector_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/math/math.hpp"
#include "vs_app_module.hpp"

// Calling program for state estimation application

/**
 * Basic constructor
 */
gridpack::voltage_stability::VSAppModule::VSAppModule(void)
{
}

/**
 * Basic destructor
 */
gridpack::voltage_stability::VSAppModule::~VSAppModule(void)
{
}

void gridpack::voltage_stability::VSAppModule::transferPFtoVS(
    boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network,
    boost::shared_ptr<VSNetwork>
    vs_network)
{
  
  int numBus = pf_network->numBuses();
  int i;
  gridpack::component::DataCollection *vsData;
  gridpack::component::DataCollection *pfData;
  double rval;
  for (i=0; i<numBus; i++) {
    pfData = pf_network->getBusData(i).get();
    vsData = vs_network->getBusData(i).get();
    pfData->getValue("BUS_PF_VMAG",&rval);
    vsData->setValue(BUS_VOLTAGE_MAG,rval);
    ///printf("Step0 bus%d mag = %f\n", i+1, rval);
    pfData->getValue("BUS_PF_VANG",&rval);
    vsData->setValue(BUS_VOLTAGE_ANG,rval);
    int ngen = 0;
    if (pfData->getValue(GENERATOR_NUMBER, &ngen)) {
      int j;
      for (j=0; j<ngen; j++) {
        pfData->getValue("GENERATOR_PF_PGEN",&rval,j);
        vsData->setValue(GENERATOR_PG,rval,j);
        pfData->getValue("GENERATOR_PF_QGEN",&rval,j);
        vsData->setValue(GENERATOR_QG,rval,j);
      }
    }
  }
}


void gridpack::voltage_stability::VSAppModule::solvepowerflow(void)
{
  p_config = gridpack::utility::Configuration::configuration();
  p_config->open(p_configfile,p_comm);

  p_configcursor = p_config->getCursor("Configuration.Powerflow");
  bool useNonLinear = false;
  useNonLinear = p_configcursor->get("UseNonLinear", useNonLinear);

  boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network(new gridpack::powerflow::PFNetwork(p_comm));
  
  p_pfapp = new gridpack::powerflow::PFAppModule;
  p_pfapp->readNetwork(pf_network, p_config);
  p_pfapp->initialize();
  if (useNonLinear) {
    p_pfapp->nl_solve();
  } else {
    p_pfapp->solve();
  }
  p_pfapp->write();
  p_pfapp->saveData();
  
  boost::shared_ptr<gridpack::voltage_stability::VSNetwork>
    vs_network(new gridpack::voltage_stability::VSNetwork(p_comm));
    
  pf_network->clone<gridpack::voltage_stability::VSBus,
		    gridpack::voltage_stability::VSBranch>(vs_network);

  transferPFtoVS(pf_network,vs_network);
}

void gridpack::voltage_stability::VSAppModule::write_FVSI()
{
  printf("\n Branch Fast Voltage Stability Index");
  printf("\n Branch Power Flow\n");
  printf("\n        Bus 1       Bus 2   CKT         P"
                  "                    Q         FVSI\n");
  gridpack::voltage_stability::VSBranch branch;
  bool solve = branch.FVSIWrite();
}
