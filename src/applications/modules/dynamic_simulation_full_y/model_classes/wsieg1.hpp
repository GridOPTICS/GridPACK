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
 * @Latested modification with control blocks: Jul 26, 2023
 * 
 * @brief  
 * 
 * 
 */

#ifndef _wsieg1_h_
#define _wsieg1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_governor_model.hpp"
#include <string>
#include "cblock.hpp"
#include "dblock.hpp"

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

    void setGV0(double gv);

    /**
     * Set internal state parameter in governor
     * @param name character string corresponding to state variable
     * @param value new value for state parameter
     * @return false if no variable corresponding to name is found
     */
    bool setState(std::string name, double value);

    /**
     * Get internal state parameter in governor
     * @param name character string corresponding to state variable
     * @param value current value for state parameter
     * @return false if no variable corresponding to name is found
     */
    bool getState(std::string name, double *value);


  private:

    // Governor WSIEG1 Parameters read from dyr
    double K, T1, T2, T3, Uo, Uc, Pmax, Pmin;
    double T4, K1, K2, T5, K3, K4, T6, K5, K6, T7, K7, K8;
    double Db1, Err, Db2;
    double Gv1, PGv1, Gv2, PGv2, Gv3, PGv3, Gv4, PGv4, Gv5, PGv5;
    double Iblock;

    double GV0;

    // Outputs: Mechnical Power Gen1 and Gen 2
    double Pmech1, Pmech2;

    bool SecondGenExists, OptionToModifyLimitsForInitialStateLimitViolation;

    double Pref;
    double w;

    Deadband Db1_blk; // is DBInt (Intentional Deadband) implemented in Deadband block?
    LeadLag Leadlag_blk;
    Integrator P_blk;
    Deadband Db2_blk; // is BackLash (Unintentional Deadband) implemented in Deadband block?
    PiecewiseSlope NGV_blk; // is GainBlock a PiecewiseSlope block? yes
    Filter Filter_blk1;
    Filter Filter_blk2;
    Filter Filter_blk3;
    Filter Filter_blk4;

    double GV;

    void computeModel(double t_inc, IntegrationStage int_flag);

};
}  // dynamic_simulation
}  // gridpack
#endif
