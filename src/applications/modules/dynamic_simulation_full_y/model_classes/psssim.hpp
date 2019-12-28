/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   psssim.hpp
 * @author Renke Huang
 * @Last modified:   Aug. 20, 2019
 * 
 * @brief  
 * 
 * 
 */

#ifndef _psssim_h_
#define _psssim_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_pss_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class PsssimModel : public BasePssModel
{
  public:
    /**
     * Basic constructor
     */
    PsssimModel();

    /**
     * Basic destructor
     */
    virtual ~PsssimModel();

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

    double getVstab();
	 
    void setOmega(double omega);
	
	double getBusFreq(int busnum);
	
	void setWideAreaFreqforPSS(double freq);	


  private:

    //PSSSIM parameters from dyr
	int inputtype;
	int bus1, bus2, bus3, bus4, bus5, bus6;
    double gaink, tw, t1, t2, t3, t4, maxout, minout;
    double kp;

    // PSSSIM state variables
    double dx1pss, dx2pss, dx3pss;
	double dx1pss_1, dx2pss_1, dx3pss_1;
	double x1pss, x2pss, x3pss;
	double x1pss_1, x2pss_1, x3pss_1;
	
	//PSS output
	double pssout_vstab;    
   
    // PSSSIM inputs
    double genspd;
	double wideareafreq;

    // Keep around 
    double psscon1, psscon2;

};
}  // dynamic_simulation
}  // gridpack
#endif
