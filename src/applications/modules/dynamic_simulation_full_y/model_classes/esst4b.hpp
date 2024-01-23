/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   esst4b.hpp
 * 
 * @brief  ESST4B
 * 
 * 
 */

#ifndef _esst4b_h_
#define _esst4b_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Esst4bModel : public BaseExciterModel
{
  public:
    /**
     * Basic constructor
     */
    Esst4bModel();

    /**
     * Basic destructor
     */
    virtual ~Esst4bModel();

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
     * FEX function
     * @ param IN
     */
    double FEX(double IN);

    /**
     * CalculateVb function
     * @ param Vterm, Theta, Ir, Ii, LadIfd
     */
    double CalculateVb(double Vterm, double Theta, double Ir, double Ii, double LadIfd);

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

    void setVuel(double vtmp);

    void setVs(double vtmp);

    /** 
     * Set the value of the terminal current
     */
    void setIri(double vIr, double vIi);

    void setVoel(double vtmp);
    
    void setVcomp(double vtmp);

  private:

    //double S10, S12; 

    // Exciter ESST4B parameters from dyr
    double Tr, Kpr, Kir, Vrmax, Vrmin, Ta, Kpm;
    double Kim, Vmmax, Vmmin, Kg, Kp, KI;
    double Vbmax, Kc, Xl, Thetap;

    bool zero_TR; // Time constant TR for measurement block is zero, no transfer function
    bool zero_TA; // Time constant TA for measurement block is zero, no transfer function  
    
    bool zero_KIM; // Integrator gain KIM for PI controller is zero, no transfer function
    bool zero_KIR; // Integrator gain KIR for PI controller is zero, no transfer function

    Filter Filter_blkR;
    PIControl PIControl_blkR;
    Filter Filter_blkA;
    PIControl PIControl_blmM;
    LVGate LVGate_blk;

    double Vuel, Vs;
    double Vref; // Reference voltage
    double Vmeas; // Ouptut of voltage measurement block
    double Voel;
    
    double Kpang, Vgmax; // TBD
   
    // ESST4B inputs
    double Vcomp, Vterm, Theta, Ir, Ii, LadIfd, Vstab; // Ir, Ii: terminal current
 
    // Field Voltage Output
    double Efd;

    // Keep around 
    double Kpvr, Kpvi, Kpir, Kpii;

    double presentMag, presentAng;
  
    bool OptionToModifyLimitsForInitialStateLimitViolation;

    void computeModel(double t_inc,IntegrationStage int_flag);

};
}  // dynamic_simulation
}  // gridpack
#endif
