/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wtara1.hpp
 * @author Shuangshuang Jin 
 * @Created on:   Nov 15, 2022
 * 
 * @brief  
 * 
 * 
 */

#ifndef _wtara1_h_
#define _wtara1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_mechanical_model.hpp"
#include <string>

namespace gridpack {
namespace dynamic_simulation {
class Wtara1Model : public BaseMechanicalModel
{
  public:
    /**
     * Basic constructor
     */
    Wtara1Model();

    /**
     * Basic destructor
     */
    virtual ~Wtara1Model();

    /**
     * Load parameters from DataCollection object into mechanical model
     * @param data collection of mechanical parameters from input files
     * @param index of mechanical on bus
     * TODO: might want to move this functionality to BaseMechanicalModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Initialize mechanical model before calculation
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
     * Set the Theta parameter inside the Aerodynamic model
     * @param Theta value of the pitch angle
     */
    void setPitchAngle(double Theta1);
    
    /**
     * Get the value of the mechanical power
     * @return value of mechanical power
     */
    double getMechanicalPower();
    /**
     * Set the waerot parameter inside the Aerodynamic model
     * @param waerot value of the turbine speed
     */
    double getWaerot();

  /**
   * Set the wt parameter inside the Aerodynamic model
   * @param wt value of the turbine speed
   */
  void  setWt(double wt1);

    
  private:

    // Mechanical WTARA1 Parameters read from dyr
    double Ka, Theta0, Pm0;
    
    // Input:
    double Theta; // External pitch angle
    double wt; // Turbine speed

    // Outputs:
    double Pmech; // Mechnical Power Gen
    double waerot; // aero-dynamic torque
};
}  // dynamic_simulation
}  // gridpack
#endif
