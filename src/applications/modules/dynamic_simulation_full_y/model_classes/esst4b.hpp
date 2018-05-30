/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   esst4b.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   Oct 12, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _esst4b_h_
#define _esst4b_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"

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

  private:

    //double S10, S12; 

    // Exciter ESST4B parameters from dyr
    double Tr, Kpr, Kir, Vrmax, Vrmin, Ta, Kpm;
    double Kim, Vmmax, Vmmin, Kg, Kp, KI;
    double Vbmax, Kc, Xl, Kpang, Vgmax;
    
    // ESST4B state variables
    double x1Vm, x2Vcomp, x3Va, x4Vr;    
    double x1Vm_1, x2Vcomp_1, x3Va_1, x4Vr_1;
    double dx1Vm, dx2Vcomp, dx3Va, dx4Vr;    
    double dx1Vm_1, dx2Vcomp_1, dx3Va_1, dx4Vr_1;
   
    // ESST4B inputs
    double Vcomp, Vterm, Theta, Ir, Ii, LadIfd, Vstab;
 
    // Field Voltage Output
    double Efd;

    // Keep around 
    double Vref;
    double Kpvr, Kpvi, Kpir, Kpii;

    double presentMag, presentAng;
  
    bool OptionToModifyLimitsForInitialStateLimitViolation;

};
}  // dynamic_simulation
}  // gridpack
#endif
