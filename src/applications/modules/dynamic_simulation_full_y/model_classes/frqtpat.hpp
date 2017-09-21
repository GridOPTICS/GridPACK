/**
 * @file   frqtpat.hpp
 * @author Renke HuANG
 * @Last modified:   July 26, 2016
 * 
 * @brief  
 * 
 * 
 */

#ifndef _frqtpat_h_
#define _frqtpat_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_relay_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class FrqtpatRelay : public BaseRelayModel {
public:
    /**
     * Basic constructor
     */
    FrqtpatRelay();

    /**
     * Basic destructor
     */
    virtual ~FrqtpatRelay();

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
     */
    void getTripStatus(int &itrip, int &itrip_prev);
	
	
  private:
	
	//parameters
	int imins, ifreq_bus;
	std::string sgenid;
	double dfreq_low, dfreq_up, dpickuptime, dbreakertime;
	
		
	// internal variables
	// internal variables
	gridpack::ComplexType pbus_volt_freq_cplx; // point to the bus voltage phasor
	double dvol_freq;		// bus voltage freq where relay is monitoring
	int icount_pickup_lowfreq;		// To count for Relay pickup time for lower frequency condition
	int icount_pickup_upfreq;		// To count for Relay pickup time for upper frequency condition
	int icount_breaker;		// To count for Breaker’s time delay
	int iflag;				//  Flag=1: counting for the breaker’s delay;
							//  Flag=0: counting for relay’s pick-up time. 
							
							
	// output variables
	int igen_trip;			//  generator_trip signal: if igen_trip == 1, take action to do load shedding.		
	int igen_trip_prev;	

};
}  // dynamic_simulation
}  // gridpac
#endif
