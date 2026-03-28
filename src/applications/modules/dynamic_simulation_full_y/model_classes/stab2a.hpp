/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   stab2a.hpp
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief  STAB2A — Simple PSS stabilizer (speed-input)
 *
 *  Signal chain:
 *    dw = omega - 1.0
 *    u1 = KT * dw                         (gain)
 *    u2 = T*s/(1+T*s) * u1                (washout)
 *    u3 = (1+T1*s)/(1+T2*s) * u2          (lead-lag 1)
 *    u4 = (1+T3*s)/(1+T4*s) * u3          (lead-lag 2)
 *    Vstab = clamp(u4, H2, H1)            (output limiter)
 *
 *  DYR format:
 *    I, 'STAB2A', ID, KT, T, T1, T2, T3, T4, H1, H2
 */

#ifndef _stab2a_h_
#define _stab2a_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_pss_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Stab2aModel : public BasePssModel
{
  public:
    Stab2aModel();
    virtual ~Stab2aModel();

    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    void init(double mag, double ang, double ts);
    void predictor(double t_inc, bool flag);
    void corrector(double t_inc, bool flag);

    double getVstab();
    void setOmega(double omega);

  private:
    // Parameters
    double KT, T, T1, T2, T3, T4, H1, H2;

    // Input
    double omega;  // rotor speed (pu)

    // Output
    double Vstab;

    // Transfer blocks
    Cblock   WO_blk;    // KT*T*s/(1+T*s)  — gain + washout combined
    LeadLag  LL1_blk;   // (1+T1*s)/(1+T2*s)
    LeadLag  LL2_blk;   // (1+T3*s)/(1+T4*s)

    // Flags
    bool zero_T2, zero_T4;

    void computeModel(double t_inc, IntegrationStage int_flag);
};
}  // dynamic_simulation
}  // gridpack
#endif
