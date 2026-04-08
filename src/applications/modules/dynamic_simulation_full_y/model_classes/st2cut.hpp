/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   st2cut.hpp
 * @author Yousu Chen
 * @Updated:  March 2026
 *
 * @brief  ST2CUT two-input PSS model (speed + power or dual-speed).
 *
 * Signal chain:
 *   sig1 = K1/(1+sT1) * input1
 *   sig2 = K2/(1+sT2) * input2
 *   IN = sig1 + sig2
 *   WO  = T3*s/(1+sT4) * IN   [washout, K=T3, T=T4]
 *   LL1 = (1+sT5)/(1+sT6) * WO
 *   LL2 = (1+sT7)/(1+sT8) * LL1
 *   LL3 = (1+sT9)/(1+sT10) * LL2
 *   Vstab = clamp(LL3, LSMIN, LSMAX)
 *
 * Input selection (MODE):
 *   1 = rotor speed deviation
 *   2 = bus frequency deviation
 *   (Only MODE=1 is currently used; input = omega - 1.0)
 */

#ifndef _st2cut_h_
#define _st2cut_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_pss_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class St2cutModel : public BasePssModel
{
  public:
    St2cutModel();
    virtual ~St2cutModel();

    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    void init(double mag, double ang, double ts);
    void predictor(double t_inc, bool flag);
    void corrector(double t_inc, bool flag);

    double getVstab();
    void setOmega(double omega);
    void setVterminal(double mag);

  private:
    // Parameters
    int    MODE, MODE2;
    double K1, K2;
    double T1, T2, T3, T4, T5, T6, T7, T8, T9, T10;
    double LSMAX, LSMIN;
    double VCU, VCL;    // voltage thresholds: PSS disabled outside [VCL, VCU]

    // Inputs
    double omega;   // rotor speed channel 1 (pu)
    double omega2;  // rotor speed channel 2 (pu)
    double Vterm;   // terminal voltage magnitude (pu)

    // Output
    double Vstab;

    // Transfer blocks
    Filter  In1_blk;  double sig1;  // K1/(1+sT1)
    Filter  In2_blk;  double sig2;  // K2/(1+sT2)
    Cblock  WO_blk;                 // T3*s/(1+sT4)
    LeadLag LL1_blk;  double LL1;  // (1+sT5)/(1+sT6)
    LeadLag LL2_blk;  double LL2;  // (1+sT7)/(1+sT8)
    LeadLag LL3_blk;  double LL3;  // (1+sT9)/(1+sT10)

    // Flags
    bool zero_T1, zero_T2, zero_T4, zero_T6, zero_T8, zero_T10;

    void computeModel(double t_inc, IntegrationStage int_flag);
};
}  // dynamic_simulation
}  // gridpack
#endif
