/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   exdc1.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   Jun 11, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _exdc1_h_
#define _exdc1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Exdc1Model : public BaseExciterModel
{
  public:
    /**
     * Basic constructor
     */
    Exdc1Model();

    /**
     * Basic destructor
     */
    virtual ~Exdc1Model();

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

    // Exciter EXDC1 parameters from dyr
    double TR, KA, TA, TB, TC, Vrmax, Vrmin;
    //double KE, TE, KF, TF1, SWITCH;
    double KE, TE, KF, TF, SWITCH; // TF?
    double E1, SE1, E2, SE2;

    // EXDC1 state variables
    double x1, x2, x3, x4, x5;
    double x1_1, x2_1, x3_1, x4_1, x5_1;
    double dx1, dx2, dx3, dx4, dx5;
    double dx1_1, dx2_1, dx3_1, dx4_1, dx5_1;

    // Field Voltage Output
    double Efd;

    // Field Current Output
    double LadIfd;

    double Vref;

    double Vterminal, w; 

    //boost::shared_ptr<BaseGeneratorModel> p_generator;
};
}  // dynamic_simulation
}  // gridpack
#endif
