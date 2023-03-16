/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   gast.hpp
 * @author Shrirang Abhyankar
 * @Last modified:   November 7, 2022
 * 
 * @brief  
 * Gas-Turbine governor model
 * 
 */

#ifndef _gast_h_
#define _gast_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_governor_model.hpp"
#include "cblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class GastModel : public BaseGovernorModel
{
  public:
    /**
     * Basic constructor
     */
    GastModel();

    /**
     * Basic destructor
     */
    virtual ~GastModel();

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

    /** 
     * Set the governor bus number
     */
    //void setExtBusNum(int ExtBusNum);
	
    /** 
     * Set the governor generator id
     */
    //void setExtGenId(std::string ExtGenId);

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

    // Governor Gast Parameters read from dyr
    double R, T1, T2, T3, AT, KT, Dt, VMAX, VMIN;

    // Gast control blocks
    Filter fuel_valve_block;
    Filter fuel_flow_block;
    Filter exh_temp_block;

    // Outputs: Mechnical Power
    double Pmech;

    // The exhaust block forms a feedback system that is hard
    // to resolve. We hold the value of exhaust block output
    // and use it in the calculation with one step delay in the
    // feedback loop
    double exh_temp_block_out;

    // Inputs
    double Pref;   // Load reference
    double delta_w;  // speed deviation

};
}  // dynamic_simulation
}  // gridpack
#endif
