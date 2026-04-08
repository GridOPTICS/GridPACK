/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   esdc2a.hpp
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief ESDC2A exciter model (IEEE DC2A with terminal-voltage-scaled limits).
 *
 * Key difference from EXDC2:
 *   VRMAX/VRMIN scale with terminal voltage Vt at every time step:
 *     VRmax_eff = (VRMAX == 0 ? 999 : VRMAX) * Vt
 *     VRmin_eff = VRMIN * Vt
 *
 * PSS/E DYR parameter order (same as EXDC2):
 *   TR  KA  TA  TB  TC  VRMAX  VRMIN  KE  TE  KF  TF1  SWITCH  E1  SE1  E2  SE2
 */

#ifndef _esdc2a_h_
#define _esdc2a_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"
#include <string>
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Esdc2aModel : public BaseExciterModel
{
  public:
    Esdc2aModel();
    virtual ~Esdc2aModel();

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
    double TR, KA, TA, TB, TC, Vrmax_param, Vrmin_param;
    double Vrmax_c;   // = (VRMAX == 0 ? 999 : VRMAX) — Vt-scaling factor
    double KE, TE, KF, TF1;
    double E1, SE1, E2, SE2;

    // Model inputs
    double Ec;    // Terminal voltage
    double Vref;  // Reference voltage
    double Vs;    // Stabilizer voltage signal

    // Model output
    double Efd;   // Field voltage

    // Model transfer blocks
    Filter     Vmeas_blk;
    double     Vmeas;

    LeadLag    Leadlag_blk;
    double     VLL;

    Filter     Regulator_blk;      // Anti-windup; limits scaled by Vt each step
    double     VR;
    GainLimiter Regulator_gain_blk;

    Cblock     Feedback_blk;
    double     VF;

    Integrator Output_blk;

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
