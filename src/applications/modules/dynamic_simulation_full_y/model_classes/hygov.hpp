/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hygov.hpp
 * @author Shrirang Abhyankar
 * @Created:   November 7, 2022
 * 
 * @brief  
 * Hydro turbine governor model
 * 
 */

#ifndef _hygov_h_
#define _hygov_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_governor_model.hpp"
#include "cblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class HygovModel : public BaseGovernorModel
{
  public:
    /**
     * Basic constructor
     */
    HygovModel();

    /**
     * Basic destructor
     */
    virtual ~HygovModel();

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

  private:

    // Governor Hygov Parameters read from dyr
    double R, r, TR, TF, TG, VELM, GMAX, GMIN, TW,AT, Dt, qNL;

    // Hygov control blocks
    Filter    filter_block; // output e
    PIControl gate_block; // output c
    Filter    opening_block; // output g
    Integrator turbine_flow_block; // output q

    double opening_block_out, gate_block_out,turbine_flow_block_out,filter_block_in;
    // Outputs: Mechnical Power
    double Pmech;

    // Inputs
    double nref;     // reference
    double delta_w;  // speed deviation

};
}  // dynamic_simulation
}  // gridpack
#endif
