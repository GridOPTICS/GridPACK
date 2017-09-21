/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   lvshbl.cpp
 * @author Renke HuANG
 * @Last modified:   July 26, 2016
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_relay_model.hpp"
#include "lvshbl.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::LvshblRelay::LvshblRelay(void)
{
	// parameters
    dloadshed_volt1 = 0.8;   // first load shedding point, p.u.
	dpickup_T1 = 0.1;		// first point pick up time, seconds
	dloadshed_frac1 = 0.0;	// first fraction of load to be shed, p.u.
	dbreakertime = 0.2;		// breaker time, seconds
	jbus_relay_load = -1;	// bus number of load bus where relay is located
	
	// internal variables
	pbus_volt_full = gridpack::ComplexType(0.0,0.0);
	dvol_mag 	   = 1.0;		// load bus voltage magnitude where relay is located
	icount_pickup  = 0;			// To count for Relay pickup time for lower voltage condition
	icount_breaker = 0;			// To count for Breaker’s time delay
	iflag = 0;						//  Flag=1: counting for the breaker’s delay;
								//  Flag=0: counting for relay’s pick-up time. 
												
	// output variables
	iload_shed = 0;			//  Load_shedding signal: if iload_shed == 1, take action to do load shedding.	
	iload_shed_prev = 0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::LvshblRelay::~LvshblRelay(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::LvshblRelay::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
	//if (!data->getValue(RELAY_IBUS, &UEL, idx)) UEL = 0.0; // 
	if (!data->getValue(RELAY_LID, &slid,           idx)) slid = ""; // TBD: LID
	if (!data->getValue(RELAY_JBUS,&jbus_relay_load,idx)) jbus_relay_load = 0; // TBD: JBUS
	if (!data->getValue(RELAY_V1, &dloadshed_volt1, idx)) dloadshed_volt1 = 0.0; // TBD: v1
	if (!data->getValue(RELAY_T1, &dpickup_T1,      idx)) dpickup_T1 = 0.0; 	// TBD: T1
	if (!data->getValue(RELAY_F1, &dloadshed_frac1, idx)) dloadshed_frac1 = 0.0; // TBD: F1
	if (!data->getValue(RELAY_V2, &dloadshed_volt2, idx)) dloadshed_volt2 = 0.0; // TBD: v2
	if (!data->getValue(RELAY_T2, &dpickup_T2,      idx)) dpickup_T2 = 0.0; 	// TBD: T2
	if (!data->getValue(RELAY_F2, &dloadshed_frac2, idx)) dloadshed_frac2 = 0.0; // TBD: F2
	if (!data->getValue(RELAY_V3, &dloadshed_volt3, idx)) dloadshed_volt3 = 0.0; // TBD: v3
	if (!data->getValue(RELAY_T3, &dpickup_T3,      idx)) dpickup_T3 = 0.0; 	// TBD: T3
	if (!data->getValue(RELAY_F3, &dloadshed_frac3, idx)) dloadshed_frac3 = 0.0; // TBD: F3
	if (!data->getValue(RELAY_TB, &dbreakertime,    idx)) dbreakertime = 0.0; // TBD: TB
        printf("dloadshed_volt1 =%8.4f \n", dloadshed_volt1);
        printf("dpickup_T1  =%8.4f \n", dpickup_T1 );
        printf("dloadshed_frac1  =%8.4f \n", dloadshed_frac1 );
        printf("dbreakertime  =%8.4f \n", dbreakertime);
	
}

/*
 *reset relay to initial status
*/
void gridpack::dynamic_simulation::LvshblRelay::ResetRelay(void)
{
	icount_pickup  = 0;			// To count for Relay pickup time for lower voltage condition
	icount_breaker = 0;			// To count for Breaker’s time delay
	iflag = 0;						//  Flag=1: counting for the breaker’s delay;
								//  Flag=0: counting for relay’s pick-up time. 
												
	// output variables
	iload_shed = 0;			//  Load_shedding signal: if iload_shed == 1, take action to do load shedding.	
	iload_shed_prev = 0;
}

/**
 * Assign to pointers to variables that need to be monitored
 * @param var vector of pointers to complex variables that are being
 * monitored
 * @return return false if not enough variable in var
 */

bool gridpack::dynamic_simulation::LvshblRelay::setMonitorVariables(std::vector<gridpack::ComplexType*> var)
{
	pbus_volt_full = *(var[0]);
	return true;
}


/**
 * Update status of relay by some time increment
 * @param delta_t time interval since last update
 * @return false if relay has been tripped
 */
bool gridpack::dynamic_simulation::LvshblRelay::updateRelay(double delta_t)
{
	///printf("\n***** Relay %d updateRelay:\n", ibus_relay_load);
	iload_shed_prev = iload_shed;
	
	// first update the voltage magnitude
	
	dvol_mag = abs(pbus_volt_full);
        printf("dvol_mag  =%8.4f \n", dvol_mag );
	
	// implement relay logic
	if (iflag == 0) {
		if (dvol_mag<dloadshed_volt1) {
			icount_pickup = icount_pickup + 1;
		}else {
			icount_pickup = 0;
		}
	}
	
	if (icount_pickup*delta_t >= dpickup_T1) {
		iflag = 1;
		icount_breaker = icount_breaker + 1;
		
		if ( icount_breaker*delta_t >= dbreakertime) {
			iload_shed = 1;
			iflag = 0;
			icount_pickup = 0;
			icount_breaker = 0;
			
		}else {
			iload_shed = 0;
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
void gridpack::dynamic_simulation::LvshblRelay::getTripStatus( int &itrip, int &itrip_prev)
{
	itrip = iload_shed;
	itrip_prev = iload_shed_prev;
}

double gridpack::dynamic_simulation::LvshblRelay::getLvshblRelayLoadFrac()
{
	return dloadshed_frac1;
}

double gridpack::dynamic_simulation::LvshblRelay::getRelayFracPar(void)
{
	return dloadshed_frac1;
}
