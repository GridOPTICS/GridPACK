/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   esst1a.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   Oct 12, 2015
 * @Last modified with control block: Sep 11, 2023
 * 
 * @brief  
 * 
 * 
 */

#ifndef _esst1a_h_
#define _esst1a_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Esst1aModel : public BaseExciterModel
{
  public:
    /**
     * Basic constructor
     */
    Esst1aModel();

    /**
     * Basic destructor
     */
    virtual ~Esst1aModel();

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

    double sqr(double x);

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
     * Set the field voltage parameter inside the exciter
     * @param fldv value of the field voltage
     */
    void setFieldVoltage(double fldv);

    /**
     * Set the field current parameter inside the exciter
     * @param fldc value of the field current
     */
    void setFieldCurrent(double fldc);

    /** 
     * Get the value of the field voltage parameter
     * @return value of field voltage
     */
    double getFieldVoltage();

    /** 
     * Get the value of the field current parameter
     * @return value of field current
     */
    double getFieldCurrent();

    /** 
     * Set the value of the Vterminal
     * @return value of field current
     */
    void setVterminal(double mag);

    /** 
     * Set the value of the Omega 
     * @return value of field current
     */
    void setOmega(double omega);

    /** 
     * Set the value of the Vcomp
     * @return value of the Vcomp
     */
    void setVcomp(double vtmp);

    void setVstab(double vtmp);

    void setVothsg(double vtmp);

    void setVuel(double vtmp);

    void setVoel(double vtmp);

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

    //double S10, S12; 

    // Exciter ESST1A parameters from dyr
    double UEL, VOS, Tr, Vimax, Vimin, Tc, Tb;
    double Tc1, Tb1, Ka, Ta, Vamax, Vamin;
    double Vrmax, Vrmin, Kc, Kf, Tf, Klr, Ilr;
   
    // ESST1A inputs
    double Vcomp, LadIfd, Vstab, Vothsg, Vuel, Voel;

    double Vterm; // Terminal voltage not Ec, Ec == Vcomp
 
    // Field Voltage Output
    double Efd;

    void computeModel(double t_inc, IntegrationStage int_flag);

    Filter Filter_blkR;
    HVGate HVGate_blk1; 
    LeadLag Leadlag_blkBC;
    LeadLag Leadlag_blkBC1;
    Filter Regulator_blk;
    GainLimiter Regulator_gain_blk;
    HVGate HVGate_blk2; 
    LVGate LVGate_blk; 
    Cblock Feedback_blkF;

    double Vf; // Output of Feedback block
    bool   zero_TA;      // Time constant TA for regulator block zero, no transfer function
    bool   zero_TR;      // Time constant TR for measurement block is zero, no transfer function

    bool zero_TF;   // Time constant TF for feedback block, if too small, no transfer function
    bool zero_TB;   // Time constant TB for first lead lag block, if too small, no transfer function
    bool zero_TB1;   // Time constant TB1 for second lead lag block, if too small, no transfer function
    bool OptionToModifyLimitsForInitialStateLimitViolation;

    double VA; // Output of Regulator blk
    double VLL1; // Output of LeadLag blk BC1
    double VLL; // Output of LeadLag blk BC
    double Vref; // Reference voltage
    double Vmeas; // Output of voltage measurement block
};
}  // dynamic_simulation
}  // gridpack
#endif
