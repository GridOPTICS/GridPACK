/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wshygp.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 11, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _wshygp_h_
#define _wshygp_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_governor_model.hpp"
#include "GainBlockClass.hpp"
#include "BackLashClass.hpp"
#include "DBIntClass.hpp"

namespace gridpack {
namespace dynamic_simulation {
class WshygpModel : public BaseGovernorModel
{
  public:
    /**
     * Basic constructor
     */
    WshygpModel();

    /**
     * Basic destructor
     */
    virtual ~WshygpModel();

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

    // Governor WSHYGP Parameters read from dyr
    double TD, KI, KD, Kp, R, Tt, KG, TP;
    double VELopen, VELclose, Pmax, Pmin, TF;
    double Trate, Aturb, Bturb, Tturb;
    
    double Db1, Err, Db2;
    double Gvx, PGvx;
    //double Iblock;

    // WSHYGP state variables
    double x1Pmech, x2Td, x3Int, x4Der, x5Pelec, x6Valve, x7Gate;
    double x1Pmech_1, x2Td_1, x3Int_1, x4Der_1, x5Pelec_1, x6Valve_1, x7Gate_1;
    double dx1Pmech, dx2Td, dx3Int, dx4Der, dx5Pelec, dx6Valve, dx7Gate;
    double dx1Pmech_1, dx2Td_1, dx3Int_1, dx4Der_1, dx5Pelec_1, dx6Valve_1, dx7Gate_1;

    // Outputs: Mechnical Power Gen1 and Gen 2
    double Pmech;

    bool SecondGenExists, OptionToModifyLimitsForInitialStateLimitViolation;

    GainBlockClass GainBlock;
    BackLashClass BackLash;
    DBIntClass DBInt;

    double Pref;
    double w;
    double GenMVABase, GenPelec;

};
}  // dynamic_simulation
}  // gridpack
#endif
