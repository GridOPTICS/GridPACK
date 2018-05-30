/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   distr1.cpp
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
#include "distr1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Distr1Relay::Distr1Relay(void)
{
    dzone1_time = 0.0; 
	dzone1_reach = 0.0; 
	dzone1_cenang = 0.0;
	dzone1_cendis = 0.0;
	dzone2_time = 0.0; 
	dzone2_reach = 0.0; 
	dzone2_cenang = 0.0;
	dzone2_cendis = 0.0;
	dsebtime = 0.0;
	dzone1_dis = 10000.0;
	dzone2_dis = 10000.0;
	
	iflag = 0;
	icount_zone1t = 0;
	icount_zone2t = 0;
	icount_breaker = 0;
	iline_trip = 0;
	iline_trip_prev = 0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Distr1Relay::~Distr1Relay(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::Distr1Relay::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
	if (!data->getValue(RELAY_ID, &sid,           idx)) sid = "";  
	if (!data->getValue(RELAY_RS, &irs,           idx)) irs = -1; 
	if (!data->getValue(RELAY_MTYPE, &imtype,     idx)) imtype = -1; 
	if (!data->getValue(RELAY_BMON, &ibmon,       idx)) ibmon = -1; 
	if (!data->getValue(RELAY_IBUS1, &ibus1,      idx)) ibus1 = -1; 
	if (!data->getValue(RELAY_JBUS1, &jbus1,      idx)) jbus1 = -1; 
	if (!data->getValue(RELAY_ID1, &sid1,         idx)) sid1 = -1; 
	if (!data->getValue(RELAY_IBUS2, &ibus2,      idx)) ibus2 = -1; 
	if (!data->getValue(RELAY_JBUS2, &jbus2,      idx)) jbus2 = -1; 
	if (!data->getValue(RELAY_ID2, &sid2,         idx)) sid2 = -1; 
	if (!data->getValue(RELAY_IBUS3, &ibus3,      idx)) ibus3 = -1; 
	if (!data->getValue(RELAY_JBUS3, &jbus3,      idx)) jbus3 = -1; 
	if (!data->getValue(RELAY_ID3, &sid3,         idx)) sid3 = -1; 
	if (!data->getValue(RELAY_ZONE1_TIME,   &dzone1_time,        idx)) dzone1_time = 0.0; 
	if (!data->getValue(RELAY_ZONE1_REACH,  &dzone1_reach,       idx)) dzone1_reach = 0.0; 
	if (!data->getValue(RELAY_ZONE1_CENANG, &dzone1_cenang,      idx)) dzone1_cenang = 0.0; 
	if (!data->getValue(RELAY_ZONE1_CENDIS, &dzone1_cendis,      idx)) dzone1_cendis = 0.0; 
	if (!data->getValue(RELAY_ZONE2_TIME,   &dzone2_time,        idx)) dzone2_time = 0.0; 
	if (!data->getValue(RELAY_ZONE2_REACH,  &dzone2_reach,       idx)) dzone2_reach = 0.0; 
	if (!data->getValue(RELAY_ZONE2_CENANG, &dzone2_cenang,      idx)) dzone2_cenang = 0.0; 
	if (!data->getValue(RELAY_ZONE2_CENDIS, &dzone2_cendis,      idx)) dzone2_cendis = 0.0; 
	if (!data->getValue(RELAY_ZONE3_TIME,   &dzone3_time,        idx)) dzone3_time = 0.0; 
	if (!data->getValue(RELAY_ZONE3_REACH,  &dzone3_reach,       idx)) dzone3_reach = 0.0; 
	if (!data->getValue(RELAY_ZONE3_CENANG, &dzone3_cenang,      idx)) dzone3_cenang = 0.0; 
	if (!data->getValue(RELAY_ZONE3_CENDIS, &dzone3_cendis,      idx)) dzone3_cendis = 0.0; 
	
	if (!data->getValue(RELAY_DIRANG, &dirang,      idx)) dirang = 0.0; 
	if (!data->getValue(RELAY_THCUR, &dthcur,       idx)) dthcur = 0.0; 
	if (!data->getValue(RELAY_SEBTIME, &dsebtime,   idx)) dsebtime = 0.0; 
	if (!data->getValue(RELAY_SERCTIME, &dserctime, idx)) dserctime = 0.0; 
	if (!data->getValue(RELAY_TRBTIME, &dtrbtime,   idx)) dtrbtime = 0.0; 
	if (!data->getValue(RELAY_TRRCTIME, &dtrrcitme, idx)) dtrrcitme = 0.0; 
	
	if (!data->getValue(RELAY_BLTYPE1, &ibltype1,    idx)) ibltype1 = -1; 
	if (!data->getValue(RELAY_BLINT1, &dblint1,      idx)) dblint1 = 0.0; 
	if (!data->getValue(RELAY_BLRO1, &dblro1,        idx)) dblro1 = 0.0; 
	if (!data->getValue(RELAY_BLTYPE2, &ibltype2,    idx)) ibltype2 = -1; 
	if (!data->getValue(RELAY_BLINT2, &dblint2,      idx)) dblint2 = 0.0; 
	if (!data->getValue(RELAY_BLRO2, &dblro2,        idx)) dblro2 = 0.0; 

	//dsebtime = 2.0; //tmp code	
	const double pi = 4.0*atan(1.0);
	c_zone1center = gridpack::ComplexType(dzone1_cendis*cos(dzone1_cenang*pi/180.0), 
										  dzone1_cendis*sin(dzone1_cenang*pi/180.0));
										  
	c_zone2center = gridpack::ComplexType(dzone2_cendis*cos(dzone2_cenang*pi/180.0), 
										  dzone2_cendis*sin(dzone2_cenang*pi/180.0));
	
	printf("dsebtime: %8.4f  \n", dsebtime);	
	printf("dserctime: %8.4f  \n", dserctime);	
	printf("dzone2_time: %8.4f  \n", dzone2_time);	
	printf("dzone2_reach: %8.4f  \n", dzone2_reach);	
	printf("dzone2_cenang: %8.4f  \n", dzone2_cenang);	
	printf("dzone2_cendis: %8.4f  \n", dzone2_cendis);	
	
	printf("c_zone1center: %8.4f + %8.4fj \n", real(c_zone1center), imag(c_zone1center));		
	printf("c_zone2center: %8.4f + %8.4fj \n", real(c_zone2center), imag(c_zone2center));			
	
}

/*
 *reset relay to initial status
*/
void gridpack::dynamic_simulation::Distr1Relay::ResetRelay(void)
{
	iflag = 0;
	icount_zone1t = 0;
	icount_zone2t = 0;
	icount_breaker = 0;
	iline_trip = 0;
	iline_trip_prev = 0;
	dzone1_dis = 10000.0;
	dzone2_dis = 10000.0;
}

/**
 * Assign to pointers to variables that need to be monitored
 * @param var vector of pointers to complex variables that are being
 * monitored
 * @return return false if not enough variable in var
 */


bool gridpack::dynamic_simulation::Distr1Relay::setMonitorVariables(std::vector<gridpack::ComplexType*> var)
{
	c_volt = *(var[0]);
	c_curr = *(var[1]);
	return true;
}

/**
 * Update status of relay by some time increment
 * @param delta_t time interval since last update
 * @return false if relay has been tripped
 */
bool gridpack::dynamic_simulation::Distr1Relay::updateRelay(double delta_t)
{
	const double pi = 4.0*atan(1.0);
	
	iline_trip_prev = iline_trip;
	
	//printRelayVoltCurr();
	
	if( iflag == 0 ) {
		dzone2_dis = abs( c_volt/c_curr - c_zone2center );
		dzone1_dis = abs( c_volt/c_curr - c_zone1center );
		
		//printf("Distr1 update: dzone1_dis = %8.4f, zone1 radio: %8.4f \n", dzone1_dis, dzone1_reach/2.0);
		//printf("Distr1 update: dzone2_dis = %8.4f, zone2 radio: %8.4f \n", dzone2_dis, dzone2_reach/2.0);
		
		if ( dzone2_dis < (dzone2_reach/2.0) ) {
			icount_zone2t++;
			//printf( "icount_zone2t: %d \n", icount_zone2t);
			
			if ( dzone1_dis < (dzone1_reach/2.0) ) {
				icount_zone1t++;
				//printf( "icount_zone1t: %d \n", icount_zone1t);
			}else {
				icount_zone1t = 0;
			}
	
		}else {
			icount_zone2t = 0;
		}	
	}
	
	if ( icount_zone1t*delta_t > dzone1_time/60 || icount_zone2t*delta_t > dzone2_time/60 ) {
		iflag = 1;
		//printf("delta_t: %8.4f \n", delta_t);
		//printf( "iflag: %d \n", iflag);
		icount_breaker++;
		//printf( "icount_breaker: %d \n", icount_breaker);
		if ( icount_breaker*delta_t > dsebtime/60.0) {
			iline_trip = 1;
			iflag = 0;
			icount_zone1t = 0;
			icount_zone2t = 0;
			icount_breaker = 0;
			
		}else {
			iline_trip = 0;
			iflag = 1;
		}
	} else {
		iflag = 0;
	}

	return true;
}

/**
 * Return internal status of relay (margin by which threshold has bee
 * exceeded)
 */
void gridpack::dynamic_simulation::Distr1Relay::getTripStatus(int &itrip, int &itrip_prev)
{
	itrip = iline_trip;
	itrip_prev = iline_trip_prev;
}

/**
 * print bus volt and current
 * 
 */
void gridpack::dynamic_simulation::Distr1Relay::printRelayVoltCurr( void)
{
	printf ("Distr1 Relay bus volt, %8.4f+%8.4fj,  bus current, %8.4f+%8.4fj,\n", real(c_volt), imag(c_volt), real(c_curr),imag(c_curr) );
}
