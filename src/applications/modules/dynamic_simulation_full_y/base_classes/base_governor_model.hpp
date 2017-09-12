/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_governor_model.hpp
 * @author Bruce Palmer
 * @Last modified:   May 18, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _base_governor_model_h_
#define _base_governor_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"

namespace gridpack {
namespace dynamic_simulation {
class BaseGovernorModel
{
  public:
    /**
     * Basic constructor
     */
    BaseGovernorModel();

    /**
     * Basic destructor
     */
    virtual ~BaseGovernorModel();

    /**
     * Load parameters from DataCollection object into governor model
     * @param data collection of governor parameters from input files
     * @param index of governor on bus
     * TODO: might want to move this functionality to BaseGovernorModel
     */
    virtual void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Initialize governor model before calculation
     * @param mag voltage magnitude
     * @param ang voltage angle
     * @param ts time step
     */
    virtual void init(double mag, double ang, double ts);

    /**
     * Predict new state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    virtual void predictor(double t_inc, bool flag);

    /**
     * Correct state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    virtual void corrector(double t_inc, bool flag);

    /**
     * Set the mechanical power parameter inside the governor
     * @param pmech value of the mechanical power
     */
    virtual void setMechanicalPower(double pmech);

    /**
     * Set the rotor speed deviation inside the governor
     * @param delta_o value of the rotor speed deviation
     */
    virtual void setRotorSpeedDeviation(double delta_o);

    /** 
     * Get the value of the mechanical power
     * @return value of mechanical power
     */
    virtual double getMechanicalPower();

    /** 
     * Get the value of the rotor speed deviation
     * @return value of rotor speed deviation
     */
    virtual double getRotorSpeedDeviation();

  private:

};
}  // dynamic_simulation
}  // gridpack
#endif
