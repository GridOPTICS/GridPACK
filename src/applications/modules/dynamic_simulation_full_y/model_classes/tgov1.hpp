/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   tgov1.hpp
 * @author Shrirang Abhyankar
 * @Last modified:   November 29, 2022
 * 
 * @brief  
 *  Turbine-governor model
 * 
 */

#ifndef _tgov1_h_
#define _tgov1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_governor_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

#define TS_THRESHOLD 4

namespace gridpack {
namespace dynamic_simulation {
class Tgov1Model : public BaseGovernorModel
{
  public:
    /**
     * Basic constructor
     */
    Tgov1Model();

    /**
     * Basic destructor
     */
    virtual ~Tgov1Model();

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
    void setRotorSpeedDeviation(double delta_w);

    /** 
     * Get the value of the mechanical power
     * @return value of mechanical power
     */
    double getMechanicalPower();

  private:

    // Governor Tgov1 Parameters read from dyr
    double R, T1, Vmax, Vmin, T2, T3, Dt;

    // Model inputs
    double Pref; // Reference power (calculated by model)
    double delta_w;     // Speed deviation (set by generator model)

    // Model output
    double Pmech;  // Mechanical power output

    // Tgov1 blocks
  LeadLag leadlag_blk;    // Lead-lag block
  double  leadlag_blk_out;  // Output of lead lag block

  Filter delay_blk; // Delay block
  double delay_blk_out; // Output of delay block
  
  void computeModel(double t_inc, IntegrationStage int_flag);

};
}  // dynamic_simulation
}  // gridpack
#endif
