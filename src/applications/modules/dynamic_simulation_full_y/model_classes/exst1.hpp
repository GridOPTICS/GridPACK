/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   exst1.hpp
 * @author Yousu Chen
 * @Updated:  March 2026
 *
 * @brief  EXST1 exciter model (PSS/E EXST1 static exciter).
 *
 * Block diagram (per PSS/E documentation):
 *   Vmeas  = 1/(1+sTR) * Vt
 *   Verr   = Vref + Vs - Vmeas - VF
 *   VI     = clamp(Verr, VIMIN, VIMAX)          [input limiter]
 *   VLL    = (1+sTC)/(1+sTB) * VI               [lead-lag; bypassed if TB=0]
 *   VR     = KA/(1+sTA) * VLL                   [anti-windup at VRMAX/VRMIN]
 *   Efd    = VR - KC * Ifd                       [no output integrator]
 *   VF     = KF * s/(1+sTF) * Efd               [rate feedback, Cblock]
 */

#ifndef _exst1_h_
#define _exst1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"
#include <string>
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Exst1Model : public BaseExciterModel
{
  public:
    Exst1Model();
    virtual ~Exst1Model();

    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    void init(double Vm, double Va, double ts);
    void predictor(double t_inc, bool flag);
    void corrector(double t_inc, bool flag);

    void setFieldVoltage(double fldv);
    double getFieldVoltage();
    void setFieldCurrent(double fldc);

    void setVterminal(double Vm);
    void setVstab(double Vstab);
    void setExtBusNum(int ExtBusNum);
    void setExtGenId(std::string ExtGenId);

    bool setState(std::string name, double value);
    bool getState(std::string name, double *value);

  private:
    // Parameters
    double TR, VIMAX, VIMIN, TC, TB, KA, TA;
    double Vrmax, Vrmin, KC, KF, TF;

    // Inputs
    double Ec;    // Terminal voltage (filtered → Vmeas)
    double Vref;  // Reference voltage (computed at init)
    double Vs;    // PSS stabilizer signal
    double LadIfd; // Field current from generator

    // Output
    double Efd;

    // Transfer blocks
    Filter   Vmeas_blk;  double Vmeas;
    LeadLag  Leadlag_blk; double VLL;
    Filter   Regulator_blk; double VR;
    GainLimiter Regulator_gain_blk;
    Cblock   Feedback_blk; double VF;

    // Flags
    bool zero_TR, zero_TB, zero_TA;

    // Bus/gen ID
    int p_bus_num;
    std::string p_gen_id;

    void computeModel(double t_inc, IntegrationStage int_flag);
};
}  // dynamic_simulation
}  // gridpack
#endif
