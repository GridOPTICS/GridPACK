/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   ieeest.hpp
 * @author Yousu Chen
 * @Updated:  March 2026
 *
 * @brief  IEEEST speed-input PSS model (IEEE 421.5-2005).
 *
 * Signal chain:
 *   dw = omega - 1.0  (speed deviation)
 *   F1a = 1/(1+sA1) * dw
 *   F1b = 1/(1+sA2) * F1a
 *   F2a = (1+sA3)/(1+sA4) * F1b
 *   F2b = (1+sA5)/(1+sA6) * F2a
 *   LL1 = (1+sT1)/(1+sT2) * F2b
 *   LL2 = (1+sT3)/(1+sT4) * LL1
 *   WO  = KS * sT5/(1+sT6) * LL2    [washout: T5 num, T6 denom]
 *   Vstab = clamp(WO, LSMIN, LSMAX)
 *
 * Note: F1 approximated as cascade of two first-order lags (standard
 * practice in time-domain simulation; exact form is 2nd-order
 * lag 1/(1+sA1+s^2*A1*A2) which differs in pole placement).
 */

#ifndef _ieeest_h_
#define _ieeest_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_pss_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class IeeeStModel : public BasePssModel
{
  public:
    IeeeStModel();
    virtual ~IeeeStModel();

    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    void init(double mag, double ang, double ts);
    void predictor(double t_inc, bool flag);
    void corrector(double t_inc, bool flag);

    double getVstab();
    void setOmega(double omega);

  private:
    // Parameters
    double A1, A2, A3, A4, A5, A6;
    double T1, T2, T3, T4, T5, T6;
    double KS, LSMAX, LSMIN;

    // Input
    double omega;  // rotor speed (pu)

    // Output
    double Vstab;

    // Transfer blocks
    Filter   F1a_blk; double F1a;  // 1/(1+sA1)
    Filter   F1b_blk; double F1b;  // 1/(1+sA2)
    LeadLag  F2a_blk; double F2a;  // (1+sA3)/(1+sA4)
    LeadLag  F2b_blk; double F2b;  // (1+sA5)/(1+sA6)
    LeadLag  LL1_blk; double LL1;  // (1+sT1)/(1+sT2)
    LeadLag  LL2_blk; double LL2;  // (1+sT3)/(1+sT4)
    Cblock   WO_blk;               // KS*sT5/(1+sT6)    [washout]

    // Flags for bypassing near-zero time constants
    bool zero_A1, zero_A2, zero_A4, zero_A6;
    bool zero_T2, zero_T4, zero_T6;

    void computeModel(double t_inc, IntegrationStage int_flag);
};
}  // dynamic_simulation
}  // gridpack
#endif
