/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   exdc2.hpp
 * @author Yousu Chen
 * @Updated:  March 2026
 *
 * @brief EXDC2 exciter model (IEEE DC2A / PSS/E EXDC2).
 *
 * PSS/E DYR parameter order (identical to EXDC1):
 *   TR  KA  TA  TB  TC  VRMAX  VRMIN  KE  TE  KF1  TF1  SWITCH  E1  SE1  E2  SE2
 */

#ifndef _exdc2_h_
#define _exdc2_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"
#include <string>
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Exdc2Model : public BaseExciterModel
{
  public:
    /**
     * Basic constructor
     */
    Exdc2Model();

    /**
     * Basic destructor
     */
    virtual ~Exdc2Model();

    /**
     * Load parameters from DataCollection object into exciter model
     * @param data collection of exciter parameters from input files
     * @param index of exciter on bus
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Saturation function
     * @param x input value
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
     */
    void setExtBusNum(int ExtBusNum);

    /**
     * Set the exciter generator id
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
    double TR, KA, TA, TB, TC, Vrmax, Vrmin;
    double KE, TE, KF1, TF1;
    double E1, SE1, E2, SE2;

    // Model inputs
    double Ec;    // Terminal voltage
    double Vref;  // Reference voltage
    double Vs;    // Stabilizer voltage signal

    // Model output
    double Efd;   // Field voltage

    // Model transfer blocks
    Filter     Vmeas_blk;          // Voltage measurement: 1/(1+sTR)
    double     Vmeas;

    LeadLag    Leadlag_blk;        // Lead-lag: (1+sTC)/(1+sTB)
    double     VLL;

    Filter     Regulator_blk;      // Anti-windup amplifier: KA/(1+sTA), x clamped to [Vrmin,Vrmax]
    double     VR;
    GainLimiter Regulator_gain_blk; // Replaces Regulator_blk when TA = 0

    Cblock     Feedback_blk;       // Washout feedback: KF1*s/(1+sTF1)
    double     VF;

    Integrator Output_blk;         // Exciter integrator: 1/(sTE)

    // Saturation curve constants
    double sat_A;
    double sat_B;
    bool   has_Sat;
    bool   has_leadlag;
    bool   zero_TA;
    bool   zero_TR;

    std::string p_gen_id;
    int p_bus_num;

    void computeModel(double t_inc, IntegrationStage int_flag);
};
}  // dynamic_simulation
}  // gridpack
#endif
