/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wtdta1_cblock.hpp
 * @author Shuangshuang Jin 
 * @Created on:   Dec 05, 2022
 * 
 * @brief  
 * 
 * 
 */

#ifndef _wtdta1_h_
#define _wtdta1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_mechanical_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"
#include <string>

namespace gridpack {
namespace dynamic_simulation {
class Wtdta1Model : public BaseMechanicalModel
{
  public:
    /**
     * Basic constructor
     */
    Wtdta1Model();

    /**
     * Basic destructor
     */
    virtual ~Wtdta1Model();

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
     * Set the Pgen parameter inside the Drive-Train model
     * @param Pgen value of the electrical power
     */
    void setPgen(double Pgen1);
    
    /**
     * Set the Pmech parameter inside the Drive-Train model
     * @param Pmech value of the mechanical power
     */
    void setPmech(double Pmech1);
    
    /**
     * Get the value of the rotational speed
     * @return value of  rotational speed
     */
    double getRotationalSpeed();

  private:

    // Mechanical Wtdta1 Parameters read from dyr:
    double Ht; // Inertia time constant
    double damp, hfrac, freq1, dshaft;
    // Initial rotational speed
    double w0;
    
    // Input: External electrical power and mechanical power
    double Pgen, Pmech;

    // Outputs: Rotational speed
    double wt;
    
    // WTDTA1 state variables
    /*double x0;
    double x1;
    double dx0;
    double dx1;*/

    Integrator itblock1, itblock2, itblock3, itblock4;
};
}  // dynamic_simulation
}  // gridpack
#endif
