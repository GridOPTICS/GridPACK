/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   vs_app.cpp
 * @author Bruce Palmer
 * @date   2018-06-20 11:07:20 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "gridpack/applications/modules/voltage_stability/vs_factory_module.hpp"
#include "vs_app_module.hpp"
//#include "vs_factory_module.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/parser/PTI33_parser.hpp"
#include "gridpack/parser/PTI34_parser.hpp"
#include "gridpack/parser/PTI35_parser.hpp"
#include "gridpack/parser/MAT_parser.hpp"
#include "gridpack/export/PSSE34Export.hpp"
#include "gridpack/export/PSSE33Export.hpp"
#include "gridpack/export/PSSE23Export.hpp"
#include "gridpack/parser/GOSS_parser.hpp"
#include "gridpack/math/math.hpp"
#include "vs_helper.hpp"
#include "gridpack/utilities/string_utils.hpp"
#include <iostream>
#include <fstream>

#define USE_REAL_VALUES

/**
 * Basic constructor
 */
gridpack::voltage_stability::VSAppModule::VSAppModule(void)
{
  current_increment = 0.0;
  p_bPVAnlyDone = true;
  zone = 0;
  sink_area = 0;
  src_area = 0;
  gt = 0.0;
  lt = 0.0;
  PV_header = false;
}

/**
 * Basic destructor
 */
gridpack::voltage_stability::VSAppModule::~VSAppModule(void)
{
}

/**
 * Increment generators real power based off specified value. 
 * Increment generators in specified area.
 * @param transfer value to increment generators real power
 * @param area index of area for incrementing generation
 * @param zone index of zone for incrementing generation
 * @param total power generation of an area
 */
void gridpack::voltage_stability::VSAppModule::IncrementGeneratorRealPower(
    double inc, int area, int zone, double gt)
{
  v_factory->IncrementGeneratorRealPower(inc,area,zone,gt);
}

/**
 * Increment load power based off specified value. 
 * Increment loads in specified area.
 * @param transfer value to increment load real power
 * @param area index of area for incrementing load
 * @param zone index of zone for incrementing load
 * @param total active power demand of the area
 */
void gridpack::voltage_stability::VSAppModule::IncrementLoadPower(
    double inc, int area, int zone, double lt)
{
  v_factory->IncrementLoadPower(inc,area,zone,lt);
}

/**
 * Return the total real power load for all loads in the zone. If zone
 * less than 1, then return the total load real power for the area
 * @param area index of area
 * @param zone index of zone
 * @return total load
 */
double gridpack::voltage_stability::VSAppModule::getTotalLoadRealPower(int area, int zone)
{
  return p_factory->getTotalLoadRealPower(area,zone);
}

/**
 * Check to see if PV Analysis is complete
 * @return return true if parameter is not found
 */
bool gridpack::voltage_stability::VSAppModule::isPVAnlyDone(double increment, double max_increment)
{
  if (current_increment >= (max_increment)){
	  p_bPVAnlyDone = true;
  }
	return p_bPVAnlyDone;
}

/**
 * Set up PV Curve internal parameters and initialize
 */
void gridpack::voltage_stability::VSAppModule::InitializePVCurve(std::string filename)
{    
  std::cout<<"1 Done "<<std::endl;
  /*
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Powerflow");
  */ 
  std::cout<<"2 Done "<<std::endl;
  //gt = v_factory->getTotalGenRealPower(src_area,zone);
  lt = getTotalLoadRealPower(sink_area, zone);
  p_bPVAnlyDone = false;    
  std::cout<<"1 Done "<<std::endl;
  /*if (!cursor->get("PVAnalysisData",&filename)) {
     printf("No PV Analysis output data file specified\n");
  }
  else{}*/    
  std::cout<<"3 Done "<<std::endl;
  static int numBus = p_network->numBuses();
  int i;
  double bus_varray[numBus] = {};
  std::ofstream file;
  file.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
  file.is_open();
  if(!PV_header){
    file << "Incremental Transfer (MW),";
    for (i=0; i<numBus; i++) {
      if (p_network->getActiveBus(i)) {
        gridpack::voltage_stability::VSBus *bus =
          dynamic_cast<gridpack::voltage_stability::VSBus*>
          (p_network->getBus(i).get());
        file << "Bus " << bus->getOriginalIndex() << ",";
      }
   }
 file << "\n";
 PV_header = true;
 file.close();
    
  }
}

/**
 * Execute one transfer increment
 */
void gridpack::voltage_stability::VSAppModule::IncrementPVCurveStep()
{
  IncrementGeneratorRealPower(increment, src_area, zone, gt);
  IncrementLoadPower(increment, sink_area, zone, lt);
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Powerflow");
  std::string filename;
  if (!cursor->get("PVAnalysisData",&filename)) {
     printf("No PV Analysis output data file specified\n");
  }
  else{
     static int numBus = p_network->numBuses();
     int i;
     double bus_varray[numBus] = {};
     std::ofstream file;
     file.open(filename.c_str(), std::ios_base::app);
     file.is_open();
     file << current_increment << ",";
     for (i=0; i<numBus; i++) {
       if (p_network->getActiveBus(i)) {
         gridpack::voltage_stability::VSBus *bus =
           dynamic_cast<gridpack::voltage_stability::VSBus*>
           (p_network->getBus(i).get());
         file << bus->getVoltage() << ",";
       }
     }
    file << "\n";
    file.close();
    
  }
  gt = gt + increment;
  lt = lt + increment;
  current_increment = current_increment + increment;
}
