/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wsieg1.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 11, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _wsieg1_h_
#define _wsieg1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_governor_model.hpp"
#include "GainBlockClass.hpp"
#include "BackLashClass.hpp"
#include "DBIntClass.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Wsieg1Model : public BaseGovernorModel
{
  public:
    /**
     * Basic constructor
     */
    Wsieg1Model();

    /**
     * Basic destructor
     */
    virtual ~Wsieg1Model();

    /**
     * Load parameters from DataCollection object into governor model
     * @param data collection of governor parameters from input files
     * @param index of governor on bus
     * TODO: might want to move this functionality to BaseGovernorModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Initialize governor model before calculation
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
     * Set the mechanical power parameter inside the governor
     * @param pmech value of the mechanical power
     */
    void setMechanicalPower(double pmech);

    /**
     * Set the rotor speed deviation inside the governor
     * @param delta_o value of the rotor speed deviation
     */
    void setRotorSpeedDeviation(double delta_o);

    /** 
     * Get the value of the mechanical power
     * @return value of mechanical power
     */
    double getMechanicalPower();

    /** 
     * Get the value of the rotor speed deviation
     * @return value of rotor speed deviation
     */
    //double getRotorSpeedDeviation();

  private:

    // Governor WSIEG1 Parameters read from dyr
    double K, T1, T2, T3, Uo, Uc, Pmax, Pmin;
    double T4, K1, K2, T5, K3, K4, T6, K5, K6, T7, K7, K8;
    double Db1, Err, Db2;
    //double Gv1, PGv1, Gv2, PGv2, Gv3, PGv3, Gv4, PGv4, Gv5, PGv5;
    double Iblock;

    // WSIEG1 state variables
    double x1LL, x2GovOut, x3Turb1, x4Turb2, x5Turb3, x6Turb4;
    double x1LL_1, x2GovOut_1, x3Turb1_1, x4Turb2_1, x5Turb3_1, x6Turb4_1;
    double dx1LL, dx2GovOut, dx3Turb1, dx4Turb2, dx5Turb3, dx6Turb4;
    double dx1LL_1, dx2GovOut_1, dx3Turb1_1, dx4Turb2_1, dx5Turb3_1, dx6Turb4_1;

    // Outputs: Mechnical Power Gen1 and Gen 2
    double Pmech1, Pmech2;

    bool SecondGenExists, OptionToModifyLimitsForInitialStateLimitViolation;

    GainBlockClass GainBlock;
    BackLashClass BackLash;
    DBIntClass DBInt;

    double Pref;
    double w;

};
}  // dynamic_simulation
}  // gridpack
#endif
