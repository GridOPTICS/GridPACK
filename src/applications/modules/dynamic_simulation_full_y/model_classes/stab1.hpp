/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   stab1.hpp
 * @author Shuangshuang Jin
 * @Last modified:   Apr 03, 2023
 * 
 * @brief  
 * 
 * 
 */

#ifndef _stab1_h_
#define _stab1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_pss_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Stab1Model : public BasePssModel
{
  public:
    /**
     * Basic constructor
     */
    Stab1Model();

    /**
     * Basic destructor
     */
    virtual ~Stab1Model();

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

    double getVothsg();
	 
    void setSpeed(double input_speed);


  private:

    //Model parameters
	  double j, j1, j2, j3, j4, j5, j6;
    double K, T, T1, T3, T2, T4, Hlim;
    //double kp;

    // Model inputs
    double speed; // Speed (pu)
    
    // Model outputs
    double vothsg;
    
    // Model parameters
    
    // Blocks
    Cblock Feedback_blk; // Feedback block
    double VF; // Output of feedback block
    
    LeadLag Leadlag_blk1, Leadlag_blk2; // Lead lag block
    double VLL1, VLL2; // Outputs of lead lag block
    
    GainLimiter Hlim_blk; // GainLimiter block
    
    /**
       computeModel - Used both by predictor and corrector
    **/
    void computeModel(double t_inc,IntegrationStage int_flag);

};
}  // dynamic_simulation
}  // gridpack
#endif
