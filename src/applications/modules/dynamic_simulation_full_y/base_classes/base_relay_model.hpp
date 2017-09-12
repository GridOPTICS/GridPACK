/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_relay_model.hpp
 * @author Bruce Palmer
 * @Last modified:   July 12, 2016
 * 
 * @brief  
 * 
 * 
 */

#ifndef _base_relay_model_h_
#define _base_relay_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"

namespace gridpack {
namespace dynamic_simulation {
class BaseRelayModel
{
  public:
    /**
     * Basic constructor
     */
    BaseRelayModel();

    /**
     * Basic destructor
     */
    virtual ~BaseRelayModel();

    /**
     * Load parameters from DataCollection object into relay model
     * @param data collection of relay parameters from input files
     * @param index of relay on bus
     */
    virtual void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Assign to pointers to variables that need to be monitored
     * @param var vector of pointers to complex variables that are being
     * monitored
     * @return return false if not enough variable in var
     */
    virtual bool setMonitorVariables(std::vector<gridpack::ComplexType*> var);

    /**
     * Update status of relay by some time increment
     * @param delta_t time interval since last update
     * @return false if relay has been tripped
     */
    virtual bool updateRelay(double delta_t);

    /**
     * Return internal status of relay (margin by which threshold has bee
     * exceeded)
     * @return current value of threshold
     */
    virtual void  getTripStatus(int &itrip, int &itrip_prev);
	bool  getOperationStatus(void);
	void  setOperationStatus( bool sta);
	
	virtual double getRelayFracPar(void);

  private:
	 bool boperationstatus;  // true: relay  included in dynamic simulation, 
							 // false: relay not included in dynamic simulation,

};
}  // dynamic_simulation
}  // gridpack
#endif
