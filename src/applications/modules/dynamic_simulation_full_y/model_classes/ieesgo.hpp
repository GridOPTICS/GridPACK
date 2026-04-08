/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ieesgo.hpp
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief  IEESGO steam turbine governor model.
 *
 */

#ifndef _ieesgo_h_
#define _ieesgo_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_governor_model.hpp"
#include <string>
#include "cblock.hpp"
#include "dblock.hpp"

#ifndef TS_THRESHOLD
#define TS_THRESHOLD 4
#endif

namespace gridpack {
namespace dynamic_simulation {
class IeesgoModel : public BaseGovernorModel
{
  public:
    IeesgoModel();
    virtual ~IeesgoModel();

    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    void init(double mag, double ang, double ts);

    void predictor(double t_inc, bool flag);
    void corrector(double t_inc, bool flag);

    void setMechanicalPower(double pmech);
    void setRotorSpeedDeviation(double delta_w);
    double getMechanicalPower();

    bool setState(std::string name, double value);
    bool getState(std::string name, double *value);

  private:
    // Parameters
    double T1, T2, T3, T4, T5, T6;
    double K1, K2, K3;
    double Pmax, Pmin;

    // I/O
    double w;       // speed deviation (set externally)
    double Pmech;   // mechanical power output
    double Pref0;   // computed reference power

    // Blocks
    Filter   F1_blk;   // K1/(1+T1*s)
    LeadLag  F2_blk;   // (1+T2*s)/(1+T3*s)
    Filter   F3_blk;   // 1/(1+T4*s)
    Filter   F4_blk;   // K2/(1+T5*s)
    Filter   F5_blk;   // K3/(1+T6*s)

    void computeModel(double t_inc, IntegrationStage int_flag);
};
}  // dynamic_simulation
}  // gridpack
#endif
