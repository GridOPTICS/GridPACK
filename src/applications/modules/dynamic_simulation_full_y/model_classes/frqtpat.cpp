/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   frqtpat.cpp
 * @author Renke HuANG
 * @Last modified:   July 26, 2016
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_relay_model.hpp"
#include "frqtpat.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::FrqtpatRelay::FrqtpatRelay(void)
{
     imins = -1; // TBD: MINS
	 ifreq_bus = -1; // TBD: FREBUS
	 sgenid = ""; // TBD: GENID
	 dfreq_low =0.0; // TBD: FL
	 dfreq_up =200.0; // TBD: FU
	 dpickuptime = 0.0; // TBD: TP
	 dbreakertime = 0.0; // TBD: TB
	 
	 pbus_volt_freq_cplx = gridpack::ComplexType(0.0,0.0);
	 dvol_freq = 60.0;
	 icount_pickup_lowfreq = 0;
	 icount_pickup_upfreq = 0;
	 icount_breaker = 0;
	 igen_trip = 0;	 
	 igen_trip_prev = 0;
     iflag = 0;	
}	

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::FrqtpatRelay::~FrqtpatRelay(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::FrqtpatRelay::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
	
	if (!data->getValue(RELAY_MINS, &imins,           idx)) imins = -1; // TBD: MINS
	if (!data->getValue(RELAY_FREBUS, &ifreq_bus,     idx)) ifreq_bus = -1; // TBD: FREBUS
	if (!data->getValue(RELAY_GENID, &sgenid,         idx)) sgenid = ""; // TBD: GENID
	if (!data->getValue(RELAY_FL, &dfreq_low,         idx)) dfreq_low =0.0; // TBD: FL
	if (!data->getValue(RELAY_FU, &dfreq_up,          idx)) dfreq_up =200.0; // TBD: FU
	if (!data->getValue(RELAY_TP, &dpickuptime,       idx)) dpickuptime = 0.0; // TBD: TP
	if (!data->getValue(RELAY_TB, &dbreakertime,      idx)) dbreakertime = 0.0; // TBD: TB
}

/*
 *reset relay to initial status
*/
void gridpack::dynamic_simulation::FrqtpatRelay::ResetRelay(void)
{
	igen_trip = 0;
	igen_trip_prev = 0;
	
	dvol_freq = 60.0;
	icount_pickup_lowfreq = 0;
	icount_pickup_upfreq = 0;
	icount_breaker = 0; 
    iflag = 0;	
	
}

/**
 * Assign to pointers to variables that need to be monitored
 * @param var vector of pointers to complex variables that are being
 * monitored
 * @return return false if not enough variable in var
 */

bool gridpack::dynamic_simulation::FrqtpatRelay::setMonitorVariables(std::vector<gridpack::ComplexType*> var)
{
	pbus_volt_freq_cplx = *(var[0]);
	return true;
}

/**
 * Update status of relay by some time increment
 * @param delta_t time interval since last update
 * @return false if relay has been tripped
 */
bool gridpack::dynamic_simulation::FrqtpatRelay::updateRelay(double delta_t)
{
	///printf("\n***** FRQTPAT Relay %d updateRelay:\n", ibus_relay_load);
	igen_trip_prev = igen_trip;
	
	// first update the voltage frequency
	
	dvol_freq = real(pbus_volt_freq_cplx);
	
	// implement relay logic
	// count for relay’s pick-up time
	if ( iflag==0 ) {
		//check under frequency condition
		if ( dvol_freq<dfreq_low ) {
			icount_pickup_lowfreq++;
		}else {
			icount_pickup_lowfreq = 0;
		}
		//check up frequency condition
		if ( dvol_freq>dfreq_up ) {
			icount_pickup_upfreq++;
		}else {
			icount_pickup_upfreq = 0;
		}
	}
	
	//check if the pick-up time is satisfied
	if ( icount_pickup_lowfreq*delta_t>=dpickuptime 
		|| icount_pickup_upfreq*delta_t>=dpickuptime ) {
		iflag = 1;
		icount_breaker++;
		
		//check if the breaker delay condition is satisfied
		if ( icount_breaker*delta_t>dbreakertime ) {
			igen_trip = 1;
			iflag = 0;
			icount_pickup_lowfreq = 0;
			icount_pickup_upfreq = 0;
			icount_breaker = 0;			
		}else {
			igen_trip = 0;
			iflag = 1;
		}	
	}else {
		iflag = 0;
	}
	
	return true;
}

/**
 * Return internal status of relay (margin by which threshold has bee
 * exceeded)
 * @return current value of threshold
 */
void gridpack::dynamic_simulation::FrqtpatRelay::getTripStatus(int &itrip, int &itrip_prev)
{
	itrip = igen_trip;
	itrip_prev = igen_trip_prev;
}
