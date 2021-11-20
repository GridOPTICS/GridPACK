/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   reeca1.hpp
 * @author Renke Huang
 * @Last modified:   May. 16, 2021
 * 
 * @brief  
 * 
 * 
 */


#ifndef _reeca1_h_
#define _reeca1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_generator_model.hpp"
#include "DelayBlock.hpp"
#include "DelayBlockwithLimit.hpp"
#include "PIBlockwithLimit.hpp"
#include "LeadLagBlock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Reeca1Model : public BaseExciterModel
{
  public:
    /**
     * Basic constructor
     */
    Reeca1Model();

    /**
     * Basic destructor
     */
    virtual ~Reeca1Model();

    /**
     * Load parameters from DataCollection object into exciter model
     * @param data collection of exciter parameters from input files
     * @param index of exciter on bus
     * TODO: might want to move this functionality to BaseExciterModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Initialize exciter model before calculation
     * @param mag voltage magnitude
     * @param ang voltage angle
     * @param ts time step 
     */
    void init(double mag, double ang, double ts);

    /**
     * Predict new state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    void predictor(double t_inc, bool flag);

    /**
     * Correct state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    void corrector(double t_inc, bool flag);

    /** 
     * Set the value of the Vterminal
     * @return value of field current
     */
    void setVterminal(double mag);
	
	void setIpcmdIqcmd(double ipcmd, double iqcmd); 
	
	void setPrefQext(double pref, double qext); 
	
	double getPref();
	double getQext();
	double getIpcmd();
	double getIqcmd();
	
	// Yuan added below 2020-6-23
	/** 
	 * Set the exciter bus number
	 * @return value of exciter bus number
	 */
	void setExtBusNum(int ExtBusNum);
	
	/** 
	 * Set the exciter generator id
	 * @return value of generator id
	 */
	void setExtGenId(std::string ExtGenId);
	
	// Yuan added above 2020-6-23
	
	// public availabe varaibles to make setting values easy
	double ipcmd, iqcmd;  
	double qext, pref;
	double Vterm, vfilt;
	
   //three state blocks
	gridpack::dynamic_simulation::DelayBlock delayblock_x1vfilt;
	gridpack::dynamic_simulation::DelayBlockwithLimit delayblocklimit_x2pord;
	gridpack::dynamic_simulation::DelayBlock delayblock_x3qv;

  private:
  
	bool OptionToModifyLimitsForInitialStateLimitViolation;

	bool bmodel_debug;
    
    double p_sbase, MVABase, trv, dbd1, dbd2, vdip, vup, kqv, imax, tiq, dpmax, dpmin, pmax, pmin, tpord;
	double vq1, iq1, vq2, iq2, vq3, iq3, vq4, iq4, vp1, ip1, vp2, ip2, vp3, ip3, vp4, ip4;
	
	double iqmax, iqmin, ipmax, ipmin;
	

    //double Vterm, Theta, Ir, Ii;
    
    double x1vfilt_0, x2pord_0, x3qv_0;
    double x1vfilt_1, x2pord_1, x3qv_1;
    double dx1vfilt_0, dx2pord_0, dx3qv_0;
    double dx1vfilt_1, dx2pord_1, dx3qv_1;

    std::string p_ckt;
    int p_bus_id;


};
}  // dynamic_simulation
}  // gridpack
#endif
