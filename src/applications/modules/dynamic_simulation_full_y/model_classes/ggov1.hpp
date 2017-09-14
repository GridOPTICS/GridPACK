/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ggov1.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 11, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _ggov1_h_
#define _ggov1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_governor_model.hpp"
#include "GainBlockClass.hpp"
#include "BackLashClass.hpp"
#include "DBIntClass.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Ggov1Model : public BaseGovernorModel
{
  public:
    /**
     * Basic constructor
     */
    Ggov1Model();

    /**
     * Basic destructor
     */
    virtual ~Ggov1Model();

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

    // Governor GGOV1 Parameters read from dyr
    double Rselect, Flag, R, Tpelec, MaxErr;
    double MinErr, Kpgov, Kigov, Kdgov, Tdgov;
    double Vmax, Vmin, Tact, Kturb, Wfnl, Tb;
    double Tc, Teng, Tfload, Kpload, Kiload;
    double Ldref, Dm, Ropen, Rclose, Kimw;
    double Aset, Ka, Ta, Trate, Db, Tsa, Tsb;
    double Rup, Rdown;
    
    double Db1, Err, Db2;
    double Gvx, PGvx;
    //double Iblock;

    // GGOV1 state variables
    double x1Pelec, x2GovDer, x3GovInt, x4Act, x5LL, x6Fload;
    double x7LoadInt, x8LoadCtrl, x9Accel, x10TempLL;
    double x1Pelec_1, x2GovDer_1, x3GovInt_1, x4Act_1, x5LL_1, x6Fload_1;
    double x7LoadInt_1, x8LoadCtrl_1, x9Accel_1, x10TempLL_1;
    double dx1Pelec, dx2GovDer, dx3GovInt, dx4Act, dx5LL, dx6Fload;
    double dx7LoadInt, dx8LoadCtrl, dx9Accel, dx10TempLL;
    double dx1Pelec_1, dx2GovDer_1, dx3GovInt_1, dx4Act_1, dx5LL_1, dx6Fload_1;
    double dx7LoadInt_1, dx8LoadCtrl_1, dx9Accel_1, dx10TempLL_1;

    // Outputs: Mechnical Power Gen1 and Gen 2
    double Pmech;

    bool SecondGenExists, OptionToModifyLimitsForInitialStateLimitViolation;

    GainBlockClass GainBlock;
    BackLashClass BackLash;
    DBIntClass DBInt;

    double Pref, Pmwset, KigovKpgov, KiLoadKpLoad, LdRefslashKturb, LastLowValueSelect, LeadLagOut;
    double w;
    double GenMVABase, GenPelec;

};
}  // dynamic_simulation
}  // gridpack
#endif
