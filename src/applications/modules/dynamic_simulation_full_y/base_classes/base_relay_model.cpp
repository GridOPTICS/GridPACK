/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_relay_model.cpp
 * @author Bruce Palmer
 * @Last modified:   July 12, 2016
 * 
 * @brief  
 * 
 * 
 */

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_relay_model.hpp"

/**
 * Basic constructor
 */
gridpack::dynamic_simulation::BaseRelayModel::BaseRelayModel()
{
	boperationstatus = true;
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::BaseRelayModel::~BaseRelayModel()
{
}

/**
 * Load parameters from DataCollection object into relay model
 * @param data collection of relay parameters from input files
 * @param index of relay on bus
 */
void gridpack::dynamic_simulation::BaseRelayModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
}

/**
 * Assign to pointers to variables that need to be monitored
 * @param var vector of pointers to complex variables that are being
 * monitored
 * @return return false if not enough variable in var
 */
bool gridpack::dynamic_simulation::BaseRelayModel::setMonitorVariables(
    std::vector<gridpack::ComplexType*> var)
{
  return true;
}

/**
 * Update status of relay by some time increment
 * @param delta_t time interval since last update
 * @return false if relay has been tripped
 */
bool gridpack::dynamic_simulation::BaseRelayModel::updateRelay(
    double delta_t)
{
  return true;
}

/**
 * Return internal status of relay (margin by which threshold has bee
 * exceeded)
 * @return current value of threshold
 */
void gridpack::dynamic_simulation::BaseRelayModel::getTripStatus(int &itrip, int &itrip_prev)
{

}

double gridpack::dynamic_simulation::BaseRelayModel::getRelayFracPar(void)
{
	return 1.0;
}

bool  gridpack::dynamic_simulation::BaseRelayModel::getOperationStatus(void)
{
	return boperationstatus;
}
void  gridpack::dynamic_simulation::BaseRelayModel::setOperationStatus( bool sta)
{
	boperationstatus = sta;
}
