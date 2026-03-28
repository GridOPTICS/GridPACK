/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ieeet2.hpp
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief  IEEET2 — IEEE Type 2 DC Exciter (PSS/E DC2A equivalent)
 *
 *  Block diagram (same as IEEET1 except two differences):
 *   1. Rate feedback: Efd → 1/(1+TF2*s) → KF*TF1*s/(1+TF1*s) → VF
 *   2. VR output limits scale with terminal voltage:
 *        VRmax_eff = VRMAX * Vt,  VRmin_eff = VRMIN * Vt
 *
 *  DYR format:
 *    I, 'IEEET2', ID, TR, KA, TA, VRMAX, VRMIN, KE, TE, KF, TF1, TF2,
 *                 E1, SE1, E2, SE2
 */

#ifndef _ieeet2_h_
#define _ieeet2_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"
#include <string>
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Ieeet2Model : public BaseExciterModel
{
  public:
    Ieeet2Model();
    virtual ~Ieeet2Model();

    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    double Sat(double x);

    void init(double Vm, double Va, double ts);
    void predictor(double t_inc, bool flag);
    void corrector(double t_inc, bool flag);

    void setFieldVoltage(double fldv);
    double getFieldVoltage();
    void setVterminal(double Vm);
    void setVstab(double vstab);
    void setExtBusNum(int ExtBusNum);
    void setExtGenId(std::string ExtGenId);
    bool setState(std::string name, double value);
    bool getState(std::string name, double *value);

  private:
    // Model parameters
    double TR, KA, TA, Vrmax, Vrmin;
    double KE, TE, KF, TF1, TF2;
    double E1, SE1, E2, SE2;

    // Model inputs
    double Ec;   // Terminal voltage
    double Vref; // Reference voltage
    double Vs;   // Stabilizer signal

    // Model output
    double Efd;

    // Transfer blocks
    Filter      Vmeas_blk;        // 1/(1+TR*s)
    double      Vmeas;

    Filter      Regulator_blk;    // KA/(1+TA*s)  with Vt-scaled limits
    double      VR;
    GainLimiter Regulator_gain_blk; // used when TA=0

    Filter      TF2_blk;          // 1/(1+TF2*s)  — inner feedback lag
    Cblock      Feedback_blk;     // KF*TF1*s/(1+TF1*s)
    double      VF;

    Integrator  Output_blk;       // 1/(TE*s)

    // Flags
    double sat_A, sat_B;
    bool   has_Sat;
    bool   zero_TA;
    bool   zero_TR;
    bool   zero_TF2;

    std::string p_gen_id;
    int p_bus_num;

    void computeModel(double t_inc, IntegrationStage int_flag);
};
}  // dynamic_simulation
}  // gridpack
#endif
