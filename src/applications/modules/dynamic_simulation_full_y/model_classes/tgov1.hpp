/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   tgov1.hpp
 * @author Renke Huang
 * @Last modified:   June 17, 20120
 * 
 * @brief  
 * 
 * 
 */

#ifndef _tgov1_h_
#define _tgov1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_governor_model.hpp"
//#include "GainBlockClass.hpp"
//#include "BackLashClass.hpp"
//#include "DBIntClass.hpp"

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

    /** 
     * Get the value of the rotor speed deviation
     * @return value of rotor speed deviation
     */
    //double getRotorSpeedDeviation();

  private:

    // Governor Tgov1 Parameters read from dyr
    double R, T1, Vmax, Vmin, T2, T3, Dt;

    // Tgov1 state variables
    double x1pow, x2val;
    double x1pow_1, x2val_1;
    double dx1pow, dx2val;
    double dx1pow_1, dx2val_1;

    // Outputs: Mechnical Power Gen1 and Gen 2
    double Pmech;

    double Pref;
    double delta_w;
	bool  bdebugmodel;

};
}  // dynamic_simulation
}  // gridpack
#endif
