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
 * 
 * @brief  
 * 
 * 
 */

#ifndef _esst1a_h_
#define _esst1a_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"

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
	
	void setVstab(double vstab);

  private:

    //double S10, S12; 

    // Exciter ESST1A parameters from dyr
    double UEL, VOS, Tr, Vimax, Vimin, Tc, Tb;
    double Tc1, Tb1, Ka, Ta, Vamax, Vamin;
    double Vrmax, Vrmin, Kc, Kf, Tf, Klr, Ilr;
    double Ta1;
    
    // ESST1A state variables
    double x1Va, x2Vcomp, x3LL1, x4LL2, x5Deriv;    
    double x1Va_1, x2Vcomp_1, x3LL1_1, x4LL2_1, x5Deriv_1;    
    double dx1Va, dx2Vcomp, dx3LL1, dx4LL2, dx5Deriv;    
    double dx1Va_1, dx2Vcomp_1, dx3LL1_1, dx4LL2_1, dx5Deriv_1;    
   
    // ESST1A inputs
    double Vcomp, Vterm, LadIfd, Vstab;
 
    // Field Voltage Output
    double Efd;

    // Keep around 
    double Vref;

    double presentMag, presentAng;
  
    bool OptionToModifyLimitsForInitialStateLimitViolation;

};
}  // dynamic_simulation
}  // gridpack
#endif
