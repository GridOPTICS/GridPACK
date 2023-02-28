/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ieeet1.hpp
 * @author Shrirang Abhyankar
 * @Updated:  December 25, 2022
 * 
 * @brief IEEET1 exciter model.
 * 
 * 
 */

#ifndef _ieeet1_h_
#define _ieeet1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"
#include <string>
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Ieeet1Model : public BaseExciterModel
{
  public:
    /**
     * Basic constructor
     */
    Ieeet1Model();

    /**
     * Basic destructor
     */
    virtual ~Ieeet1Model();

    /**
     * Load parameters from DataCollection object into exciter model
     * @param data collection of exciter parameters from input files
     * @param index of exciter on bus
     * TODO: might want to move this functionality to BaseExciterModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Saturation function
     * @ param x
     */
    double Sat(double x);

    /**
     * Initialize exciter model before calculation
     * @param Vm voltage magnitude
     * @param Va voltage angle
     * @param ts time step 
     */
    void init(double Vm, double Va, double ts);

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
     * Set the field voltage parameter inside the exciter
     * Only used during initialization
     * @param fldv value of the field voltage
     */
    void setFieldVoltage(double fldv);

    /** 
     * Get the value of the field voltage parameter
     * @return value of field voltage
     */
    double getFieldVoltage();

    /** 
     * Set terminal voltage
     */
    void setVterminal(double Vm);

  /**
   * Set stabilizer signal
   */
    void setVstab(double vstab);
	
    /** 
     * Set the exciter bus number
     * @return value of exciter bus number
    **/
    void setExtBusNum(int ExtBusNum);
	
    /** 
     * Set the exciter generator id
     * 
    */
    void setExtGenId(std::string ExtGenId);

    /** 
     * Set internal state parameter in exciter
     * @param name character string corresponding to state variable
     * @param value new value for state parameter
     * @return false if no variable corresponding to name is found
     */
    bool setState(std::string name, double value);

    /** 
     * Get internal state parameter in exciter
     * @param name character string corresponding to state variable
     * @param value current value for state parameter
     * @return false if no variable corresponding to name is found
     */
    bool getState(std::string name, double *value);

  private:

    // Model parameters
    double TR, KA, TA, Vrmax, Vrmin;
    double KE, TE, KF, TF1;
    double E1, SE1, E2, SE2;

    // Model inputs
    double Ec; // Terminal voltage
    double Vref; // Reference voltage
    double Vs; // Stabilizer voltage signal

    // Model output
    double Efd;  // Field voltage

    // Model transfer blocks
    Filter Vmeas_blk; // Voltage measurement block
    double Vmeas; // Output of voltage measurement block
  
    Filter Regulator_blk; // Voltage regulator block
    double VR;  // Voltage regulator output
    GainLimiter Regulator_gain_blk; // Replaces Regulator block if TA = 0
  
    Cblock Feedback_blk; // Feedback block
    double VF; // Output of feedback block

    Integrator Output_blk; // Ouput block

    // Variables
  // A and B parameters for quadratic saturation curve
  double sat_A;
  double sat_B; 
  bool   has_Sat;     // Saturation active?
  bool   zero_TA;      // Time constant TA for regulator block zero, no transfer function
  bool   zero_TR;      // Time constant TR for measurement block is zero, no transfer function
  
    std::string p_gen_id; // id of the generator where the exciter is installed on
    int p_bus_num;  // bus number of the generator

  void computeModel(double t_inc,IntegrationStage int_flag);


};
}  // dynamic_simulation
}  // gridpack
#endif
