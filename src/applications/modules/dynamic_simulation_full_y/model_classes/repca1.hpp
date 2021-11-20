/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   repca1.hpp
 * @author Renke Huang
 * @Last modified:   Aug. 20, 2019
 * 
 * @brief  
 * 
 * 
 */

#ifndef _repca1_h_
#define _repca1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_plant_model.hpp"
#include "DelayBlock.hpp"
#include "DelayBlockwithLimit.hpp"
#include "LeadLagBlock.hpp"
#include "PIBlockwithLimit.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Repca1Model : public BasePlantControllerModel
{
  public:
    /**
     * Basic constructor
     */
    Repca1Model();

    /**
     * Basic destructor
     */
    virtual ~Repca1Model();

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

    void setGenPQV(double P, double Q, double Vt);

	void setBusFreq(double freq);
	
	void setPrefQext(double Pref, double Qext);
	
	void setExtBusNum(int ExtBusNum);
	
	void setExtGenId(std::string ExtGenId);
	
	double getPref( );
	
	double getQext( );
	
	double pref, qext, vref, freqref, plant_pref, genP, genQ, Vterm, busfreq;
	
	//blocks------------------------------------

    gridpack::dynamic_simulation::DelayBlock delayblock_x1vmeas;
	//gridpack::dynamic_simulation::DelayBlock delayblock_x2qmeas;
	gridpack::dynamic_simulation::PIBlockwithLimit piblock_x3qpi;
	gridpack::dynamic_simulation::LeadLagBlock leadlagblock_x4qext;
	gridpack::dynamic_simulation::DelayBlock delayblock_x5pmeas;
	gridpack::dynamic_simulation::PIBlockwithLimit piblock_x6ppi;
	gridpack::dynamic_simulation::DelayBlock delayblock_x7pref;


  private:

    bool bmodel_debug;
	
	double kc, tfltr, dbd1, dbd2, emax, emin, qmax, qmin, kp, ki, tft, tfv, fdbd1, fdbd2, ddn, dup, tp, femax, femin;
	double kpg, kig, pmax, pmin, tg;
	
	std::string p_ckt;
    int p_bus_id;

};
}  // dynamic_simulation
}  // gridpack
#endif
