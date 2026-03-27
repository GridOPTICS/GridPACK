/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   esst3a.hpp
 *
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief  ESST3A static exciter model (IEEE Std 421.5-2016 Type ST3A).
 *
 * Signal flow:
 *   Vmeas  = 1/(1+sTR) * Vt
 *   Verr   = Vref + Vs - Vmeas
 *   VI     = clamp(Verr, VIMIN, VIMAX)
 *   VLL    = (1+sTC)/(1+sTB) * VI            [lead-lag; bypassed if TB=0]
 *   VR     = KA/(1+sTA) * VLL                [anti-windup VRMIN/VRMAX]
 *   VG     = clamp(KG * Efd, 0, VGMAX)
 *   vrs    = VR - VG
 *   VM     = KM/(1+sTM) * vrs                [anti-windup VMMIN/VMMAX]
 *   VE     = |KPC*(Vd+jVq) + j*(KI+KPC*XL)*(Id+jIq)|  where KPC=KP*exp(j*THETAP)
 *   IN     = KC*LadIfd/VE
 *   FEX    = piecewise commutation factor
 *   VB     = clamp(VE*FEX, 0, VBMAX)
 *   Efd    = VM * VB
 */

#ifndef _esst3a_h_
#define _esst3a_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"
#include <string>
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Esst3aModel : public BaseExciterModel
{
  public:
    Esst3aModel();
    virtual ~Esst3aModel();

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

    void setVdqIdq(double Vd, double Vq, double Id, double Iq);

    bool setState(std::string name, double value);
    bool getState(std::string name, double *value);

  private:
    // Parameters
    double TR, VIMAX, VIMIN, KM, TC, TB, KA, TA;
    double Vrmax, Vrmin, KG, KP, KI, VBMAX, KC, XL, VGMAX, THETAP, TM;
    double Vmmax, Vmmin;

    // Derived
    double THETAP_rad;  // THETAP converted from degrees to radians
    double KPC_re, KPC_im;  // KPC = KP * exp(j*THETAP_rad) in rectangular form
    double KI_eff;           // KI + KPC_re * XL (imaginary component of VE formula)

    // Inputs
    double Ec;      // Terminal voltage magnitude
    double Vref;    // Reference voltage (computed at init)
    double Vs;      // PSS stabilizer signal
    double LadIfd;  // Field current from generator
    double Vd, Vq, Id, Iq;  // d-q terminal voltage and current

    // Output
    double Efd;

    // Internal signals
    double VR;    // Regulator output
    double VM;    // Inner field regulator output
    double VB;    // Rectifier output
    double Vmeas; // Filtered terminal voltage
    double VLL;   // Lead-lag output

    // Transfer blocks
    Filter      Vmeas_blk;
    LeadLag     Leadlag_blk;
    Filter      Regulator_blk;
    GainLimiter Regulator_gain_blk;
    Filter      InnerReg_blk;
    GainLimiter InnerReg_gain_blk;

    // Flags
    bool zero_TR, zero_TB, zero_TA, zero_TM;

    // Bus/gen ID
    int p_bus_num;
    std::string p_gen_id;

    // Compute FEX piecewise commutation factor
    double computeFex(double IN) const;

    void computeModel(double t_inc, IntegrationStage int_flag);
};
}  // dynamic_simulation
}  // gridpack
#endif
