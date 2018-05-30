/**
 * @file   lvshbl.hpp
 * @author Renke HuANG
 * @Last modified:   July 26, 2016
 * 
 * @brief  
 * 
 * 
 */

#ifndef _lvshbl_h_
#define _lvshbl_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_relay_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class LvshblRelay : public BaseRelayModel {
public:
    /**
     * Basic constructor
     */
    LvshblRelay();

    /**
     * Basic destructor
     */
    virtual ~LvshblRelay();

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
		
	double getLvshblRelayLoadFrac ();
	
	double getRelayFracPar(void);

  private:
	
	//parameters
	double dloadshed_volt1;   // first load shedding point, p.u.
	double dpickup_T1;		// first point pick up time, seconds
	double dloadshed_frac1;	// first fraction of load to be shed, p.u.
	
	double dloadshed_volt2;   // second load shedding point, p.u.
	double dpickup_T2;		// second point pick up time, seconds
	double dloadshed_frac2;	// second fraction of load to be shed, p.u.
	
	double dloadshed_volt3;   // third load shedding point, p.u.
	double dpickup_T3;		// third point pick up time, seconds
	double dloadshed_frac3;	// third fraction of load to be shed, p.u.
	
	double dbreakertime;		// breaker time, seconds
	int jbus_relay_load;	    // bus number of load bus where relay is located
	
	std::string slid;
	
	// internal variables
	gridpack::ComplexType pbus_volt_full; // point to the bus voltage phasor
	double dvol_mag;		// load bus voltage magnitude where relay is located
	int icount_pickup;		// To count for Relay pickup time for lower voltage condition
	int icount_breaker;		// To count for Breaker’s time delay
	int iflag;				//  Flag=1: counting for the breaker’s delay;
							//  Flag=0: counting for relay’s pick-up time. 
							
							
	// output variables
	int iload_shed;			//  Load_shedding signal: if iload_shed == 1, take action to do load shedding.		
	int iload_shed_prev;	//  Load_shedding signal of the previous time steps.			


};
}  // dynamic_simulation
}  // gridpac
#endif
