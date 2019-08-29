/**
 * @file   frqtpat.hpp
 * @author Renke HuANG
 * @Last modified:   July 26, 2016
 * 
 * @brief  
 * 
 * 
 */

#ifndef _distr1_h_
#define _distr1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_relay_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Distr1Relay : public BaseRelayModel {
public:
    /**
     * Basic constructor
     */
    Distr1Relay();

    /**
     * Basic destructor
     */
    virtual ~Distr1Relay();

    /**
     * Load parameters from DataCollection object into relay model
     * @param data collection of relay parameters from input files
     * @param index of relay on bus
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);
		
	/*
	 * reset relay to initial status
	 */
	void ResetRelay(void);	

    /**
     * Assign to pointers to variables that need to be monitored
     * @param var vector of pointers to complex variables that are being
     * monitored
     * @return return false if not enough variable in var
     */
    bool setMonitorVariables(std::vector<gridpack::ComplexType*> var);

    /**
     * Update status of relay by some time increment
     * @param delta_t time interval since last update
     * @return false if relay has been tripped
     */
    bool updateRelay(double delta_t);

    /**
     * Return internal status of relay (margin by which threshold has bee
     * exceeded)
     * @return current value of threshold
     */
    void getTripStatus(int &itrip, int &itrip_prev);
	
	/**
     * print bus volt and current
	 * 
	 */
	void printRelayVoltCurr (void);
	
	
  private:
	
	//parameters
	std::string sid, sid1, sid2, sid3;
	int irs, imtype, ibmon, ibus1, jbus1, ibus2, jbus2, ibus3, jbus3, ibltype1, ibltype2;
	double dzone1_time, dzone1_reach, dzone1_cenang, dzone1_cendis;
	double dzone2_time, dzone2_reach, dzone2_cenang, dzone2_cendis;
	double dzone3_time, dzone3_reach, dzone3_cenang, dzone3_cendis;
	double dirang, dthcur, dsebtime, dserctime, dtrbtime, dtrrcitme;
	double dblint1, dblro1, dblint2, dblro2;
	
	
	// internal variables
	int iflag;				//  Flag=1: counting for the breaker’s delay;
							//  Flag=0: counting for relay’s pick-up time. 
	gridpack::ComplexType c_volt, c_curr;	// voltage and current monitored by the relay		
	int icount_zone1t, icount_zone2t; //	To count time inside of Zone 1 and zone 2
	int icount_breaker;		// To count for Breaker’s time delay
	gridpack::ComplexType c_zone1center, c_zone2center; // zone 1 and zone 2 center
	double dzone1_dis, dzone2_dis;
							
	// output variables
	int iline_trip;
	int iline_trip_prev;
					


};
}  // dynamic_simulation
}  // gridpac
#endif
